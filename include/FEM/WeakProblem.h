/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * Distributed under the GNU General Public License (GPL) v2.
 */
#ifndef model_WeakProblem_H_
#define model_WeakProblem_H_
#include <WeakProblemSolvers.h>
#include <Eigen/SVD>
#include <LowRankBilinear.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace model
{
    struct WeakProblemOptions
    {
        bool useDirectSolver=true;
        double tolerance=1.e-10;
        double constraintTolerance=1.e-10;
        Eigen::Index maxIterations=1000;
    };

    struct WeakProblemResult
    {
        double relativeResidual=0.0;
        Eigen::Index iterations=0;
        bool reusedFactorization=false;
    };

    // A snapshot of the bilinear coefficients is reused until reassembleMatrix().
    // RHS values and constraints are read afresh at every solve. Not thread safe.
    // Forms with lowRankPenaltyRows() use an indefinite bordered system; select
    // solvers supporting indefinite matrices (e.g. SparseLU and BiCGSTAB, not LLT).
    template<class DirectSolverType,class IterativeSolverType>
    class WeakProblem
    {
    public:
        using Matrix=typename DirectSolverType::MatrixType;
        using Vector=Eigen::VectorXd;
        using Triplet=Eigen::Triplet<double>;
        static_assert(std::is_same_v<Matrix,typename IterativeSolverType::MatrixType>,
                      "WeakProblem solvers must use the same matrix type.");
        static_assert(std::is_same_v<typename Matrix::Scalar,double>);

    private:
        struct Constraints
        {
            std::vector<Eigen::Index> reduced;
            Vector offset;
            Eigen::Index size=0;
            Eigen::MatrixXd modes,affineNullspace,globalRows;
            Vector affineOffset,globalValues;
            std::map<size_t,double> dirichlet;
        };
        const WeakProblemOptions options;
        std::shared_ptr<const void> fieldIdentity;
        std::function<std::vector<Triplet>()> assemble;
        std::function<Vector()> force;
        std::function<Eigen::MatrixXd()> penalty;
        std::function<Constraints()> constraints;
        std::function<void(const Vector&)> commit;
        std::vector<Triplet> fullTriplets;
        std::vector<Eigen::Index> mapping;
        Eigen::MatrixXd cachedModes;
        Matrix reducedA;
        Matrix solverA;
        Eigen::MatrixXd fullPenalty;
        Eigen::MatrixXd reducedPenalty;
        Vector reducedB;
        Eigen::MatrixXd affineModes;
        Vector affineSolution;
        Eigen::MatrixXd linearRows;
        Eigen::MatrixXd reducedRows;
        Vector linearValues;
        DirectSolverType direct;
        IterativeSolverType iterative;
        bool assembled=false;
        bool factorized=false;
        size_t factorizations=0;

        static bool sameMatrix(const Eigen::MatrixXd& a,const Eigen::MatrixXd& b)
        {
            return a.rows()==b.rows() && a.cols()==b.cols() && (a.array()==b.array()).all();
        }

        template<class Field>
        Constraints makeConstraints(const Field& field) const
        {
            const auto& fe=field.fe();
            constexpr int components=Field::dofPerNode;
            const Eigen::Index n=field.gSize();
            const Eigen::Index k=affineModes.cols();
            if(k>0 && affineModes.rows()!=n)
            {
                throw std::invalid_argument("Affine FEM modes have the wrong number of DOFs.");
            }
            Constraints result;
            result.offset=Vector::Zero(n);
            result.modes=Eigen::MatrixXd::Zero(n,k);
            if(k>0)
            {
                result.modes=affineModes;
            }
            result.reduced.resize(n,-1);
            std::vector<size_t> root(n);
            for(size_t node=0;node<fe.nodeSize();++node)
            {
                const size_t representative=fe.periodicRepresentative(node);
                if(representative>=fe.nodeSize())
                {
                    throw std::runtime_error("Invalid FEM periodic representative.");
                }
                const auto jump=(field.periodicGradient()*fe.periodicShift(node)).eval();
                for(int c=0;c<components;++c)
                {
                    root[node*components+c]=representative*components+c;
                    result.offset(node*components+c)=jump(c);
                    if(field.periodicNodalOffset().size())
                    {
                        result.offset(node*components+c)+=field.periodicNodalOffset()(node*components+c)
                            -field.periodicNodalOffset()(representative*components+c);
                    }
                }
            }
            const Vector prescribedOffset=result.offset;
            std::map<size_t,size_t> anchors;
            for(const auto& bc : field.dirichletConditions())
            {
                if(bc.first>=root.size() || !std::isfinite(bc.second))
                {
                    throw std::invalid_argument("Invalid FEM Dirichlet DOF or value.");
                }
                result.dirichlet.emplace(bc.first,bc.second);
                anchors.emplace(root[bc.first],bc.first);
            }
            Eigen::MatrixXd coupling=Eigen::MatrixXd::Zero(result.dirichlet.size()-anchors.size(),k);
            Vector values=Vector::Zero(coupling.rows());
            Eigen::Index row=0;
            for(const auto& bc : result.dirichlet)
            {
                const size_t anchor=anchors.at(root[bc.first]);
                if(anchor!=bc.first)
                {
                    // Two prescribed copies in one periodic class constrain the
                    // affine coefficients, rather than conflicting automatically.
                    if(k>0)
                    {
                        coupling.row(row)=affineModes.row(bc.first)-affineModes.row(anchor);
                    }
                    values(row)=bc.second-prescribedOffset(bc.first)-result.dirichlet.at(anchor)+prescribedOffset(anchor);
                    ++row;
                }
            }
            result.affineOffset=Vector::Zero(k);
            result.affineNullspace=Eigen::MatrixXd::Identity(k,k);
            if(k>0 && coupling.rows()>0)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(coupling,Eigen::ComputeThinU|Eigen::ComputeFullV);
                svd.setThreshold(options.constraintTolerance);
                result.affineOffset=svd.solve(values);
                result.affineNullspace=svd.matrixV().rightCols(k-svd.rank());
            }
            const Vector mismatch=coupling*result.affineOffset-values;
            for(Eigen::Index i=0;i<mismatch.size();++i)
            {
                const double scale=std::max({1.0,std::abs(values(i)),coupling.row(i).norm()*result.affineOffset.norm()});
                if(std::abs(mismatch(i))>options.constraintTolerance*scale)
                {
                    throw std::invalid_argument("Incompatible affine/periodic and nodal Dirichlet conditions.");
                }
            }
            std::map<size_t,Eigen::Index> independent;
            for(Eigen::Index i=0;i<n;++i)
            {
                const auto anchor=anchors.find(root[i]);
                if(anchor!=anchors.end())
                {
                    const size_t j=anchor->second;
                    result.offset(i)+=result.dirichlet.at(j)-prescribedOffset(j);
                    if(k>0)
                    {
                        result.modes.row(i)-=affineModes.row(j);
                    }
                }
                else
                {
                    const auto entry=independent.emplace(root[i],result.size);
                    if(entry.second)
                    {
                        ++result.size;
                    }
                    result.reduced[i]=entry.first->second;
                }
            }
            result.offset+=result.modes*result.affineOffset;
            result.modes=(result.modes*result.affineNullspace).eval();
            result.size+=result.modes.cols();
            if(!result.offset.allFinite() || !result.modes.allFinite())
            {
                throw std::invalid_argument("Nonfinite FEM constraint elimination.");
            }

            // Check global affine-loading constraints against the eliminated
            // Dirichlet conditions before factorization. Remove consistent
            // redundant rows, otherwise their multipliers make the system singular.
            result.globalRows=Eigen::MatrixXd::Zero(0,result.size);
            result.globalValues=Vector::Zero(0);
            if(linearRows.rows()>0)
            {
                if(linearRows.cols()!=n || linearValues.size()!=linearRows.rows())
                {
                    throw std::invalid_argument("Invalid global FEM constraint size.");
                }
                Eigen::MatrixXd rows=Eigen::MatrixXd::Zero(linearRows.rows(),result.size);
                for(Eigen::Index i=0;i<n;++i)
                {
                    if(result.reduced[i]>=0)
                    {
                        rows.col(result.reduced[i])+=linearRows.col(i);
                    }
                }
                if(result.modes.cols()>0)
                {
                    rows.rightCols(result.modes.cols())=linearRows*result.modes;
                }
                const Vector originalOffset=linearRows*result.offset;
                Vector rhs=linearValues-originalOffset;
                Vector normalizedValues=linearValues;
                Vector normalizedOffset=originalOffset;
                for(Eigen::Index i=0;i<rows.rows();++i)
                {
                    const double scale=linearRows.row(i).norm();
                    if(scale>0.0)
                    {
                        rows.row(i)/=scale;
                        rhs(i)/=scale;
                        normalizedValues(i)/=scale;
                        normalizedOffset(i)/=scale;
                    }
                }
                Eigen::MatrixXd left=Eigen::MatrixXd::Identity(rows.rows(),rows.rows());
                Eigen::Index rank=0;
                if(result.size>0)
                {
                    Eigen::JacobiSVD<Eigen::MatrixXd> svd(rows,Eigen::ComputeFullU);
                    const double largest=svd.singularValues().size()>0 ? svd.singularValues()(0) : 0.0;
                    if(largest>0.0)
                    {
                        // An eliminated constraint can leave roundoff-sized rows;
                        // compare against the original normalized row scale as well.
                        svd.setThreshold(options.constraintTolerance*std::max(1.0,rows.norm())/largest);
                        rank=svd.rank();
                    }
                    left=svd.matrixU();
                }
                if((left.rightCols(left.cols()-rank).transpose()*rhs).norm()>
                   options.constraintTolerance*std::max({1.0,normalizedValues.norm(),normalizedOffset.norm()}))
                {
                    throw std::invalid_argument("Incompatible global affine constraints and nodal Dirichlet conditions.");
                }
                result.globalRows=left.leftCols(rank).transpose()*rows;
                result.globalValues=left.leftCols(rank).transpose()*rhs;
            }
            return result;
        }

        bool prepare(const Constraints& c)
        {
            if(!assembled)
            {
                auto entries=assemble();
                for(const auto& entry : entries)
                {
                    if(entry.row()<0 || entry.col()<0 || entry.row()>=c.offset.size() ||
                       entry.col()>=c.offset.size() || !std::isfinite(entry.value()))
                    {
                        throw std::invalid_argument("Invalid FEM matrix triplet.");
                    }
                }
                fullTriplets=std::move(entries);
                assembled=true;
            }
            if(factorized && mapping==c.reduced && sameMatrix(cachedModes,c.modes) && sameMatrix(reducedRows,c.globalRows))
            {
                return true;
            }
            factorized=false;
            fullPenalty=penalty();
            if(fullPenalty.cols()!=c.offset.size() || !fullPenalty.allFinite())
            {
                throw std::invalid_argument("Invalid low-rank FEM penalty.");
            }
            reducedPenalty=Eigen::MatrixXd::Zero(fullPenalty.rows(),c.size);
            for(Eigen::Index i=0;i<c.offset.size();++i)
            {
                if(c.reduced[i]>=0)
                {
                    reducedPenalty.col(c.reduced[i])+=fullPenalty.col(i);
                }
            }
            const Eigen::Index modes=c.modes.cols();
            const Eigen::Index nodal=c.size-modes;
            if(modes>0)
            {
                reducedPenalty.rightCols(modes)=fullPenalty*c.modes;
            }
            if(!reducedPenalty.allFinite())
            {
                throw std::runtime_error("Nonfinite reduced FEM penalty.");
            }
            std::vector<Triplet> entries;
            entries.reserve(fullTriplets.size());
            for(const auto& entry : fullTriplets)
            {
                const auto row=c.reduced[entry.row()];
                const auto col=c.reduced[entry.col()];
                if(row>=0 && col>=0)
                {
                    entries.emplace_back(row,col,entry.value());
                }
            }
            if(modes>0)
            {
                Eigen::MatrixXd action=Eigen::MatrixXd::Zero(c.offset.size(),modes);
                Eigen::MatrixXd transposeAction=action;
                for(const auto& entry : fullTriplets)
                {
                    action.row(entry.row())+=entry.value()*c.modes.row(entry.col());
                    transposeAction.row(entry.col())+=entry.value()*c.modes.row(entry.row());
                }
                for(Eigen::Index i=0;i<c.offset.size();++i)
                {
                    if(c.reduced[i]>=0)
                    {
                        for(Eigen::Index j=0;j<modes;++j)
                        {
                            entries.emplace_back(c.reduced[i],nodal+j,action(i,j));
                            entries.emplace_back(nodal+j,c.reduced[i],transposeAction(i,j));
                        }
                    }
                }
                const Eigen::MatrixXd block=c.modes.transpose()*action;
                for(Eigen::Index i=0;i<modes;++i)
                {
                    for(Eigen::Index j=0;j<modes;++j)
                    {
                        entries.emplace_back(nodal+i,nodal+j,block(i,j));
                    }
                }
            }
            reducedA.resize(c.size,c.size);
            reducedA.setFromTriplets(entries.begin(),entries.end());
            reducedA.prune([](Eigen::Index,Eigen::Index,double value)
            {
                return value!=0.0;
            });
            reducedA.makeCompressed();
            for(Eigen::Index k=0;k<reducedA.nonZeros();++k)
            {
                if(!std::isfinite(reducedA.valuePtr()[k]))
                {
                    throw std::runtime_error("Nonfinite reduced FEM matrix.");
                }
            }
            // [K  B^T; B -I] eliminates to K+B^T B without forming a dense outer product.
            entries.clear();
            for(Eigen::Index outer=0;outer<reducedA.outerSize();++outer)
            {
                for(typename Matrix::InnerIterator it(reducedA,outer);it;++it)
                {
                    entries.emplace_back(it.row(),it.col(),it.value());
                }
            }
            for(Eigen::Index row=0;row<reducedPenalty.rows();++row)
            {
                entries.emplace_back(c.size+row,c.size+row,-1.0);
                for(Eigen::Index col=0;col<c.size;++col)
                {
                    if(reducedPenalty(row,col)!=0.0)
                    {
                        entries.emplace_back(c.size+row,col,reducedPenalty(row,col));
                        entries.emplace_back(col,c.size+row,reducedPenalty(row,col));
                    }
                }
            }
            reducedRows=c.globalRows;
            for(Eigen::Index i=0;i<reducedRows.rows();++i)
            {
                for(Eigen::Index j=0;j<c.size;++j)
                {
                    entries.emplace_back(c.size+reducedPenalty.rows()+i,j,reducedRows(i,j));
                    entries.emplace_back(j,c.size+reducedPenalty.rows()+i,reducedRows(i,j));
                }
            }
            const Eigen::Index systemSize=c.size+reducedPenalty.rows()+reducedRows.rows();
            solverA.resize(systemSize,systemSize);
            solverA.setFromTriplets(entries.begin(),entries.end());
            solverA.makeCompressed();
            if(c.size>0)
            {
                if(solverA.nonZeros()==0)
                {
                    throw std::runtime_error("WeakProblem has a zero matrix and unconstrained DOFs.");
                }
                if(options.useDirectSolver)
                {
                    direct.compute(solverA);
                    if(direct.info()!=Eigen::Success)
                    {
                        throw std::runtime_error("WeakProblem direct factorization failed; check constraints and null modes.");
                    }
                }
                else
                {
                    iterative.compute(solverA);
                    if(iterative.info()!=Eigen::Success)
                    {
                        throw std::runtime_error("WeakProblem iterative setup failed.");
                    }
                }
                ++factorizations;
            }
            mapping=c.reduced;
            cachedModes=c.modes;
            factorized=true;
            return false;
        }

    public:
        template<class BWF,class LWF>
        WeakProblem(std::shared_ptr<BWF> bwf,std::shared_ptr<LWF> lwf,
                    WeakProblemOptions opts={}) : options(opts)
        {
            if(!bwf || !(options.tolerance>0.0) || !std::isfinite(options.tolerance) ||
               !(options.constraintTolerance>=0.0) || !std::isfinite(options.constraintTolerance) ||
               options.maxIterations<=0)
            {
                throw std::invalid_argument("Invalid WeakProblem form or options.");
            }
            fieldIdentity=bwf->field().ownership();
            assemble=[bwf]() { return bwf->globalTriplets(); };
            penalty=[bwf]() { return lowRankBilinearRows(*bwf); };
            constraints=[this,bwf]() { return makeConstraints(bwf->field()); };
            commit=[bwf](const Vector& x) { bwf->field().setFullDofVector(x); };
            iterative.setTolerance(options.tolerance);
            iterative.setMaxIterations(options.maxIterations);
            setRightHandSide(std::move(lwf));
        }
        WeakProblem(const WeakProblem&)=delete;
        WeakProblem(WeakProblem&&)=delete;
        WeakProblem& operator=(const WeakProblem&)=delete;
        WeakProblem& operator=(WeakProblem&&)=delete;

        template<class LWF>
        void setRightHandSide(std::shared_ptr<LWF> lwf)
        {
            if(!lwf || lwf->field().ownership()!=fieldIdentity)
            {
                throw std::invalid_argument("WeakProblem forms must refer to the same trial instance.");
            }
            force=[lwf]() { return lwf->globalVector(); };
        }

        // Additional displacement modes are independent unknowns appended to the masters.
        // Keep them fixed across RHS changes to reuse the factorization.
        void setAffineModes(const Eigen::MatrixXd& modes)
        {
            if(!modes.allFinite())
            {
                throw std::invalid_argument("Nonfinite affine FEM modes.");
            }
            affineModes=modes;
            factorized=false;
        }
        void setLinearConstraints(const Eigen::MatrixXd& rows)
        {
            if(!rows.allFinite())
            {
                throw std::invalid_argument("Nonfinite global FEM constraints.");
            }
            linearRows=rows;
            linearValues=Vector::Zero(rows.rows());
            factorized=false;
        }
        void setLinearConstraintValues(const Vector& values)
        {
            if(values.size()!=linearRows.rows() || !values.allFinite())
            {
                throw std::invalid_argument("Invalid global FEM constraint values.");
            }
            linearValues=values;
        }
        const Vector& affineCoefficients() const
        {
            return affineSolution;
        }

        // Invalidate the snapshot after changing a bilinear coefficient.
        // The next solve reassembles and recomputes the factorization.
        void reassembleMatrix()
        {
            assembled=false;
            factorized=false;
        }

        WeakProblemResult solveWithMasterSlaveConstraints()
        {
            const auto c=constraints();
            const Vector b=force();
            if(b.size()!=c.offset.size() || !b.allFinite())
            {
                throw std::invalid_argument("Invalid FEM right-hand side.");
            }
            WeakProblemResult result;
            result.reusedFactorization=prepare(c) && c.size>0;
            Vector lifted=b;
            for(const auto& entry : fullTriplets)
            {
                lifted(entry.row())-=entry.value()*c.offset(entry.col());
            }
            reducedB=Vector::Zero(c.size);
            for(Eigen::Index i=0;i<b.size();++i)
            {
                if(c.reduced[i]>=0)
                {
                    reducedB(c.reduced[i])+=lifted(i);
                }
            }
            if(c.modes.cols()>0)
            {
                reducedB.tail(c.modes.cols())=c.modes.transpose()*lifted;
            }
            if(!reducedB.allFinite())
            {
                throw std::runtime_error("Nonfinite reduced FEM force vector.");
            }
            const Vector penaltyOffset=fullPenalty*c.offset;
            Vector solverB(c.size+fullPenalty.rows()+reducedRows.rows());
            solverB.head(c.size)=reducedB;
            solverB.segment(c.size,fullPenalty.rows())=-penaltyOffset;
            if(reducedRows.rows()>0)
            {
                solverB.tail(reducedRows.rows())=c.globalValues;
            }
            reducedB-=reducedPenalty.transpose()*penaltyOffset;
            if(!solverB.allFinite() || !reducedB.allFinite())
            {
                throw std::runtime_error("Nonfinite low-rank FEM force correction.");
            }
            Vector solution=Vector::Zero(solverB.size());
            if(c.size>0)
            {
                if(options.useDirectSolver)
                {
                    solution=direct.solve(solverB);
                    if(direct.info()!=Eigen::Success)
                    {
                        throw std::runtime_error("WeakProblem direct solve failed.");
                    }
                }
                else
                {
                    solution=iterative.solve(solverB);
                    result.iterations=iterative.iterations();
                    if(iterative.info()!=Eigen::Success)
                    {
                        throw std::runtime_error("WeakProblem iterative solve did not converge.");
                    }
                }
                result.relativeResidual=(solverA*solution-solverB).norm()/std::max(1.0,solverB.norm());
                if(!solution.allFinite() || !std::isfinite(result.relativeResidual) ||
                   result.relativeResidual>options.tolerance)
                {
                    throw std::runtime_error("WeakProblem residual check failed.");
                }
            }
            Vector full=c.offset;
            for(Eigen::Index i=0;i<full.size();++i)
            {
                if(c.reduced[i]>=0)
                {
                    full(i)+=solution(c.reduced[i]);
                }
            }
            const Vector coefficients=solution.segment(c.size-c.modes.cols(),c.modes.cols());
            if(c.modes.cols()>0)
            {
                full+=c.modes*coefficients;
            }
            if(linearRows.rows()>0 && (linearRows*full-linearValues).norm()>options.tolerance*std::max(1.0,linearValues.norm()))
            {
                throw std::runtime_error("Global FEM constraint residual check failed.");
            }
            for(const auto& bc : c.dirichlet)
            {
                if(std::abs(full(bc.first)-bc.second)>options.constraintTolerance*std::max(1.0,std::abs(bc.second)))
                {
                    throw std::runtime_error("Nodal Dirichlet residual check failed.");
                }
            }
            const Vector recovered=c.affineOffset+c.affineNullspace*coefficients;
            commit(full); // No field mutation until assembly and solve have succeeded.
            affineSolution=recovered;
            return result;
        }

        // Diagnostic only: materializing a global low-rank penalty can be dense.
        Matrix reducedMatrix() const
        {
            if(reducedPenalty.rows()==0)
            {
                return reducedA;
            }
            const Eigen::MatrixXd product=reducedPenalty.transpose()*reducedPenalty;
            return Matrix(reducedA+Matrix(product.sparseView()));
        }
        const Matrix& solverMatrix() const
        {
            return solverA;
        }
        const Vector& reducedForceVector() const
        {
            return reducedB;
        }
        size_t factorizationCount() const
        {
            return factorizations;
        }
    };

    template<class DirectSolverType=DefaultWeakProblemDirectSolver<Eigen::SparseMatrix<double>>,
             class IterativeSolverType=Eigen::BiCGSTAB<Eigen::SparseMatrix<double>>,
             class BWF,class LWF>
    auto weakProblem(std::shared_ptr<BWF> bwf,std::shared_ptr<LWF> lwf,WeakProblemOptions options={})
    {
        return std::make_shared<WeakProblem<DirectSolverType,IterativeSolverType>>(std::move(bwf),std::move(lwf),options);
    }
}
#endif
