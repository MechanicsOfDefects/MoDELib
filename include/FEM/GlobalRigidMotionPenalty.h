/* This file is part of MoDELib, distributed under the GNU GPL v2. */
#ifndef model_GlobalRigidMotionPenalty_H_
#define model_GlobalRigidMotionPenalty_H_
#include <BilinearWeakExpression.h>
#include <Eigen/Eigenvalues>
#include <Quadrature.h>
#include <GaussLegendre.h>
#include <BarycentricTraits.h>
#include <memory>
#include <Eigen/SparseCore>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace model
{
    // Integral projection onto homogeneous rigid modes admitted by the constraints.
    // Only six global moments are stored; WeakProblem handles F^T F as a low-rank term.
    template<class Trial>
    class GlobalRigidMotionPenalty : public BilinearWeakExpression<GlobalRigidMotionPenalty<Trial>>
    {
    public:
        using TrialFunctionType=Trial;
        using Element=typename Trial::ElementType;
        using RigidMatrix=Eigen::Matrix<double,3,6>;
        using GramMatrix=Eigen::Matrix<double,6,6>;
        using QuadratureType=Quadrature<3,4,GaussLegendre>;
        static_assert(Trial::dim==3 && Trial::dofPerNode==3);
    private:
        const std::shared_ptr<Trial> trial;
        const double factor;
        const bool freeAffineJumps;
        Eigen::Vector3d centroid=Eigen::Vector3d::Zero();
        double length=0.0;
        double volume=0.0;
        Eigen::MatrixXd moments;
        GramMatrix gram=GramMatrix::Zero();

        static Eigen::Matrix3d crossMatrix(const Eigen::Vector3d& x)
        {
            Eigen::Matrix3d result;
            result<<0,-x(2),x(1),x(2),0,-x(0),-x(1),x(0),0;
            return result;
        }
        RigidMatrix rigid(const Eigen::Vector3d& x) const
        {
            RigidMatrix result;
            result.leftCols<3>().setIdentity();
            result.rightCols<3>()=-crossMatrix((x-centroid)/length);
            return result;
        }
        Eigen::Matrix<double,6,Trial::dofPerElement> momentKernel(const Eigen::Vector3d& a,const Element& e) const
        {
            const auto bary=BarycentricTraits<3>::x2l(a);
            return rigid(e.position(bary)).transpose()*trial->sfm(e,bary)*e.absJ(bary)/volume;
        }
        GramMatrix gramKernel(const Eigen::Vector3d& a,const Element& e) const
        {
            const auto bary=BarycentricTraits<3>::x2l(a);
            const auto R=rigid(e.position(bary));
            return R.transpose()*R*e.absJ(bary)/volume;
        }
        Eigen::MatrixXd admissibleParameters() const
        {
            GramMatrix constraints=GramMatrix::Zero();
            const auto& pairs=trial->fe().periodicNodePairs();
            if(!freeAffineJumps)
            {
                for(const auto& pair : pairs)
                {
                    // A rigid rotation omega is periodic iff omega cross p == 0.
                    RigidMatrix jump=RigidMatrix::Zero();
                    jump.rightCols<3>()=-crossMatrix(pair.shift.normalized());
                    constraints+=jump.transpose()*jump/double(pairs.size());
                }
            }
            const auto& conditions=trial->dirichletConditions();
            for(const auto& entry : conditions)
            {
                if(entry.first>=trial->gSize())
                {
                    throw std::invalid_argument("Invalid rigid-motion Dirichlet DOF.");
                }
                const Eigen::Matrix<double,1,6> row=rigid(trial->fe().node(entry.first/3).P0).row(entry.first%3);
                constraints+=row.transpose()*row/(row.squaredNorm()*double(conditions.size()));
            }
            Eigen::SelfAdjointEigenSolver<GramMatrix> eig(constraints);
            if(eig.info()!=Eigen::Success)
            {
                throw std::runtime_error("Rigid-motion constraint decomposition failed.");
            }
            const double tolerance=1.e-10*std::max(1.0,eig.eigenvalues().maxCoeff());
            const Eigen::Index count=(eig.eigenvalues().array()<tolerance).count();
            return eig.eigenvectors().leftCols(count);
        }
    public:
        GlobalRigidMotionPenalty(std::shared_ptr<Trial> u,double alpha,bool freeAffineJumpsIn=false) : trial(std::move(u)),factor(alpha),freeAffineJumps(freeAffineJumpsIn)
        {
            if(!trial || !std::isfinite(factor) || factor<=0.0)
            {
                throw std::invalid_argument("Rigid-motion penalty factor must be positive.");
            }
            const auto& fe=trial->fe();
            for(const auto& pair : fe.elements())
            {
                volume+=pair.second.simplex.vol0;
                centroid+=pair.second.simplex.vol0*pair.second.simplex.center();
            }
            if(!(volume>0.0))
            {
                throw std::runtime_error("Empty rigid-motion penalty domain.");
            }
            centroid/=volume;
            length=(fe.mesh.xMax()-fe.mesh.xMin()).norm();
            moments=Eigen::MatrixXd::Zero(6,trial->gSize());
            for(const auto& pair : fe.elements())
            {
                const auto& e=pair.second;
                Eigen::Matrix<double,6,Trial::dofPerElement> local=Eigen::Matrix<double,6,Trial::dofPerElement>::Zero();
                QuadratureType::integrate(this,local,&GlobalRigidMotionPenalty::momentKernel,e);
                QuadratureType::integrate(this,gram,&GlobalRigidMotionPenalty::gramKernel,e);
                for(int node=0;node<Trial::nodesPerElement;++node)
                {
                    moments.middleCols(3*e.node(node).gID,3)+=local.middleCols(3*node,3);
                }
            }
        }
        const auto& field() const
        {
            return trial->field();
        }
        // Sparse part only. The global term is provided by lowRankPenaltyRows().
        std::vector<Eigen::Triplet<double>> globalTriplets() const
        {
            return {};
        }
        Eigen::MatrixXd lowRankPenaltyRows() const
        {
            const Eigen::MatrixXd S=admissibleParameters();
            if(S.cols()==0)
            {
                return Eigen::MatrixXd(0,moments.cols());
            }
            const Eigen::MatrixXd G=S.transpose()*gram*S;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(G);
            if(eig.info()!=Eigen::Success || eig.eigenvalues().minCoeff()<=0.0)
            {
                throw std::runtime_error("Degenerate rigid-motion Gram matrix.");
            }
            const Eigen::MatrixXd normalization=eig.eigenvectors()*eig.eigenvalues().array().sqrt().inverse().matrix().asDiagonal()*eig.eigenvectors().transpose();
            return std::sqrt(factor)*normalization*S.transpose()*moments;
        }
        Eigen::MatrixXd rigidModes() const
        {
            const Eigen::MatrixXd S=admissibleParameters();
            Eigen::MatrixXd result(trial->gSize(),S.cols());
            for(size_t node=0;node<trial->fe().nodeSize();++node)
            {
                result.middleRows(3*node,3)=rigid(trial->fe().node(node).P0)*S;
            }
            return result;
        }
    };
}
#endif
