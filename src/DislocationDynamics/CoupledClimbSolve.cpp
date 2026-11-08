#include <CoupledClimbSolve.h>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <Utilities/SolverStageTimer.h>
namespace model
{
    void CoupledClimbSolve::compute(const Eigen::MatrixXd& K,const Eigen::SparseMatrix<double>& local,
                                    const Eigen::SparseMatrix<double>& Z)
    {
        previousSolution.resize(0);
        ready=false;
        physicalOperator={};
        directFallback.reset();
        physicalSize=K.rows();
        if(K.cols()!=physicalSize || Z.rows()!=physicalSize || local.rows()!=physicalSize
           || local.cols()!=physicalSize || !K.allFinite())
        {
            throw std::runtime_error("Invalid coupled climb matrix dimensions or entries.");
        }
        SolverStageTimer setup("Climb periodic matrix reduction ("+std::to_string(physicalSize)+" -> "+std::to_string(Z.cols())+")");
        const Eigen::MatrixXd reduced=Z.transpose()*(K*Z);
        setup.start("Climb active-row selection");
        std::vector<Eigen::Triplet<double>> selection;
        std::vector<Eigen::Triplet<double>> inactive;
        for(Eigen::Index i=0;i<reduced.rows();++i)
        {
            if(reduced.row(i).squaredNorm()>0.0 || reduced.col(i).squaredNorm()>0.0)
            {
                selection.emplace_back(i,selection.size(),1.0);
            }
            else
            {
                inactive.emplace_back(i,inactive.size(),1.0);
            }
        }
        Eigen::SparseMatrix<double> inactiveSelection(Z.cols(),inactive.size());
        inactiveSelection.setFromTriplets(inactive.begin(),inactive.end());
        inactiveMap=Z*inactiveSelection;
        Eigen::SparseMatrix<double> S(Z.cols(),selection.size());
        S.setFromTriplets(selection.begin(),selection.end());
        setup.start("Climb active matrix reduction and scaling");
        map=Z*S;
        matrix=S.transpose()*(reduced*S);
        scale.resize(matrix.rows());
        for(Eigen::Index i=0;i<matrix.rows();++i)
        {
            const double d=std::abs(matrix(i,i));
            if(!(d>0.0))
            {
                throw std::runtime_error("Active climb degree of freedom has zero self interaction.");
            }
            scale(i)=1.0/std::sqrt(d);
        }
        matrix=scale.asDiagonal()*matrix*scale.asDiagonal();
        if(matrix.rows())
        {
            setup.start("Climb local preconditioner assembly");
            Eigen::SparseMatrix<double> P=map.transpose()*local*map;
            P=scale.asDiagonal()*P*scale.asDiagonal();
            setup.start("Climb sparse preconditioner factorization (n="+std::to_string(matrix.rows())+")");
            localSolver.compute(P);
            if(localSolver.info()!=Eigen::Success)
            {
                throw std::runtime_error("Local climb preconditioner factorization failed.");
            }
        }
        ready=true;
    }

    Eigen::VectorXd CoupledClimbSolve::product(const Eigen::VectorXd& x) const
    {
        if(!physicalOperator)return matrix*x;
        const Eigen::VectorXd y=physicalOperator(map*(scale.asDiagonal()*x));
        if(y.size()!=physicalSize||!y.allFinite())throw std::runtime_error("Invalid matrix-free climb product");
        return scale.asDiagonal()*(map.transpose()*y);
    }
    void CoupledClimbSolve::computeOperator(std::function<Eigen::VectorXd(const Eigen::VectorXd&)> op,
        const Eigen::SparseMatrix<double>& local,const Eigen::SparseMatrix<double>& Z)
    {
        previousSolution.resize(0);
        ready=false;directFallback.reset();matrix.resize(0,0);physicalOperator=std::move(op);physicalSize=local.rows();
        if(!physicalOperator||local.cols()!=physicalSize||Z.rows()!=physicalSize)throw std::invalid_argument("Invalid matrix-free climb dimensions");
        SolverStageTimer reduction("Climb setup: constraint reduction and scaling",true);
        const Eigen::SparseMatrix<double> reduced=Z.transpose()*local*Z;
        std::vector<bool> active(reduced.rows(),false);
        for(int j=0;j<reduced.outerSize();++j)for(Eigen::SparseMatrix<double>::InnerIterator it(reduced,j);it;++it)
        {if(!std::isfinite(it.value()))throw std::invalid_argument("Nonfinite climb preconditioner");if(it.value()!=0){active[it.row()]=true;active[it.col()]=true;}}
        std::vector<Eigen::Triplet<double>> yes,no;
        for(Eigen::Index i=0;i<reduced.rows();++i){if(active[i])yes.emplace_back(i,yes.size(),1);else no.emplace_back(i,no.size(),1);}
        Eigen::SparseMatrix<double> S(Z.cols(),yes.size()),I(Z.cols(),no.size());S.setFromTriplets(yes.begin(),yes.end());I.setFromTriplets(no.begin(),no.end());map=Z*S;inactiveMap=Z*I;
        Eigen::SparseMatrix<double> P=S.transpose()*reduced*S;scale.resize(P.rows());
        for(Eigen::Index i=0;i<P.rows();++i){const double d=std::abs(P.coeff(i,i));if(!(d>0))throw std::runtime_error("Active FMM climb row has zero local self interaction");scale(i)=1/std::sqrt(d);}
        P=scale.asDiagonal()*P*scale.asDiagonal();
        reduction.finish();
        if(P.rows())
        {
            SolverStageTimer factorization("Climb setup: preconditioner factorization (n="+std::to_string(P.rows())+")",true);
            localSolver.compute(P);
            if(localSolver.info()!=Eigen::Success)
            {
                throw std::runtime_error("FMM climb sparse preconditioner factorization failed");
            }
        }
        ready=true;
    }

    Eigen::VectorXd CoupledClimbSolve::solve(const Eigen::VectorXd& force) const
    {
        if(!ready || force.size()!=physicalSize || !force.allFinite())
        {
            throw std::runtime_error("Climb solve requires a current finite interaction matrix and force.");
        }
        if((inactiveMap.transpose()*force).norm()>1.e-12*force.norm())
        {
            throw std::runtime_error("Nonzero climb drive on an inactive interaction row.");
        }
        const Eigen::VectorXd projected=map.transpose()*force;
        // Inactive rows are constrained to zero. The caller supplies no physical
        // drive there; all active reduced equations are checked below.
        const Eigen::VectorXd rhs=scale.asDiagonal()*projected;
        const double rhsNorm=rhs.stableNorm();
        const double projectedNorm=projected.stableNorm();
        if(rhsNorm==0.0)
        {
            previousSolution=Eigen::VectorXd::Zero(scale.size());
            return Eigen::VectorXd::Zero(physicalSize);
        }
        // Exclusive wall-clock buckets. Include warm-start and final residual
        // products, not just Arnoldi products. The matrix-free callback fences
        // device work before returning, so operator time includes completion.
        const bool profileEnabled=profilingEnabled();
        using ProfileClock=std::chrono::steady_clock;
        const auto profileStart=profileEnabled?ProfileClock::now():ProfileClock::time_point{};
        double operatorSeconds=0,arnoldiOperatorSeconds=0,orthogonalizationSeconds=0,preconditionerSeconds=0;
        std::size_t operatorCalls=0,arnoldiCalls=0,preconditionerCalls=0;
        const auto timedProduct=[&](const Eigen::VectorXd& value,bool arnoldi=false)->Eigen::VectorXd
        {
            const auto start=profileEnabled?ProfileClock::now():ProfileClock::time_point{};
            Eigen::VectorXd result=product(value);
            const double elapsed=profileEnabled?std::chrono::duration<double>(ProfileClock::now()-start).count():0.0;
            operatorSeconds+=elapsed;
            ++operatorCalls;
            if(arnoldi) { ++arnoldiCalls; arnoldiOperatorSeconds+=elapsed; }
            return result;
        };
        const Eigen::Index n=scale.size();
        const int restart=std::min<Eigen::Index>(50,n);
        const int maxIterations=std::max<Eigen::Index>(500,2*n);
        const double tolerance=1.e-10;
        const double physicalTolerance=1.e-8;
        // ||S^-1 r|| <= ||r|| / min(S). Use this bound for the Arnoldi
        // estimate, but verify both actual residuals before accepting a solve.
        const double estimatedTarget=std::min(tolerance*rhsNorm,
                                             physicalTolerance*scale.minCoeff()*projectedNorm);
        Eigen::VectorXd x=Eigen::VectorXd::Zero(n);
        Eigen::VectorXd initialResidual=rhs;
        bool warmStart=false;
        if(!directFallback && previousSolution.size()==n && previousSolution.allFinite())
        {
            const Eigen::VectorXd candidateResidual=rhs-timedProduct(previousSolution);
            // A radically different or tiny RHS can make the previous solution
            // worse than zero. Compare actual residuals in both scalings.
            if(candidateResidual.allFinite() && candidateResidual.stableNorm()<=rhsNorm &&
               candidateResidual.cwiseQuotient(scale).stableNorm()<=projectedNorm)
            {
                x=previousSolution;
                initialResidual=candidateResidual;
                warmStart=true;
            }
        }
        Eigen::VectorXd cachedResidual=initialResidual;
        bool residualCurrent=true;
        const auto actualResidual=[&]() -> const Eigen::VectorXd&
        {
            if(!residualCurrent)
            {
                cachedResidual=rhs-timedProduct(x);
                residualCurrent=true;
            }
            return cachedResidual;
        };
        int iterations=0;
        const auto failure=[&](const char* reason)
        {
            const Eigen::VectorXd& r=actualResidual();
            std::ostringstream message;
            message<<std::scientific<<std::setprecision(6)
                   <<reason<<" iterations="<<iterations<<"/"<<maxIterations
                   <<", reduced size="<<n
                   <<", scaled relative residual="<<r.stableNorm()/rhsNorm
                   <<" (target "<<tolerance<<")"
                   <<", physical relative residual="<<r.cwiseQuotient(scale).stableNorm()/projectedNorm
                   <<" (target "<<physicalTolerance<<")";
            return std::runtime_error(message.str());
        };
        if(directFallback)
        {
            SolverStageTimer cached("Climb cached LU substitution",true);
            x=directFallback->solve(rhs);
            residualCurrent=false;
        }
        // Restarted right-preconditioned GMRES, with reorthogonalization.
        if(!directFallback)
        {
            if(profilingEnabled())
        {
            std::cout<<"\n  GMRES n="<<n<<", restart="<<restart<<", max iterations="<<maxIterations
                     <<", scaled tolerance="<<tolerance<<", physical tolerance="<<physicalTolerance
                     <<", initial guess="<<(warmStart?"previous solution":"zero")
                     <<", initial scaled relative residual="<<initialResidual.stableNorm()/rhsNorm
                     <<", Eigen configured threads="<<Eigen::nbThreads()<<std::flush;
        }
            SolverStageTimer gmres("Climb GMRES",true);
            auto lastProgress=std::chrono::steady_clock::now();
            while(!directFallback && iterations<maxIterations)
            {
                const Eigen::VectorXd& residual=actualResidual();
                const double beta=residual.stableNorm();
                if(!residual.allFinite())
                {
                    throw failure("Nonfinite coupled climb residual.");
                }
                if(beta<=tolerance*rhsNorm &&
                   residual.cwiseQuotient(scale).stableNorm()<=physicalTolerance*projectedNorm)
                {
                    break;
                }
                Eigen::MatrixXd V=Eigen::MatrixXd::Zero(n,restart+1);
                Eigen::MatrixXd W=Eigen::MatrixXd::Zero(n,restart);
                Eigen::MatrixXd H=Eigen::MatrixXd::Zero(restart+1,restart);
                Eigen::VectorXd g=Eigen::VectorXd::Zero(restart+1);
                Eigen::VectorXd cs=Eigen::VectorXd::Zero(restart);
                Eigen::VectorXd sn=Eigen::VectorXd::Zero(restart);
                V.col(0)=residual/beta;
                g(0)=beta;
                int count=0;
                for(int j=0;j<restart && iterations<maxIterations;++j)
                {
                    ++iterations;
                    const auto preconditionerStart=profileEnabled?ProfileClock::now():ProfileClock::time_point{};
                    W.col(j)=localSolver.solve(V.col(j));
                    if(profileEnabled)
                    {
                        preconditionerSeconds+=std::chrono::duration<double>(ProfileClock::now()-preconditionerStart).count();
                    }
                    ++preconditionerCalls;
                    Eigen::VectorXd w=timedProduct(W.col(j),true);
                    const auto orthogonalizationStart=profileEnabled?ProfileClock::now():ProfileClock::time_point{};
                    for(int pass=0;pass<2;++pass)
                    {
                        for(int i=0;i<=j;++i)
                        {
                            const double h=V.col(i).dot(w);
                            H(i,j)+=h;
                            w-=h*V.col(i);
                        }
                    }
                    H(j+1,j)=w.norm();
                    const bool breakdown=H(j+1,j)<=1.e-14;
                    if(!breakdown)
                    {
                        V.col(j+1)=w/H(j+1,j);
                    }
                    if(profileEnabled)
                    {
                        orthogonalizationSeconds+=std::chrono::duration<double>(ProfileClock::now()-orthogonalizationStart).count();
                    }
                    for(int i=0;i<j;++i)
                    {
                        const double h=cs(i)*H(i,j)+sn(i)*H(i+1,j);
                        H(i+1,j)=-sn(i)*H(i,j)+cs(i)*H(i+1,j);
                        H(i,j)=h;
                    }
                    const double norm=std::hypot(H(j,j),H(j+1,j));
                    if(!(norm>0.0))
                    {
                        throw failure("Coupled climb GMRES breakdown.");
                    }
                    cs(j)=H(j,j)/norm;
                    sn(j)=H(j+1,j)/norm;
                    H(j,j)=norm;
                    H(j+1,j)=0.0;
                    g(j+1)=-sn(j)*g(j);
                    g(j)=cs(j)*g(j);
                    count=j+1;
                    const auto now=std::chrono::steady_clock::now();
                    if(std::chrono::duration<double>(now-lastProgress).count()>=5.0)
                    {
                        std::ostringstream progress;
                        progress<<"\n    GMRES: "<<iterations<<" iterations, estimated relative residual "
                                <<std::scientific<<std::setprecision(2)<<std::abs(g(j+1))/rhsNorm;
                        std::cout<<progress.str()<<std::endl;
                        lastProgress=now;
                    }
                    if(std::abs(g(j+1))<=estimatedTarget || breakdown)
                    {
                        break;
                    }
                }
                x+=W.leftCols(count)*H.topLeftCorner(count,count).triangularView<Eigen::Upper>().solve(g.head(count));
                residualCurrent=false;
            }
            const Eigen::VectorXd& finalResidual=actualResidual();
            if(profilingEnabled())
        {
            std::cout<<"\n  GMRES finished: iterations="<<iterations
                     <<", scaled relative residual="<<finalResidual.stableNorm()/rhsNorm
                     <<", physical relative residual="<<finalResidual.cwiseQuotient(scale).stableNorm()/projectedNorm<<std::flush;
        }
        else
        {
            std::ostringstream summary;
            summary<<"\n    GMRES: "<<iterations<<" iterations, final relative residual "
                   <<std::scientific<<std::setprecision(2)<<finalResidual.stableNorm()/rhsNorm;
            std::cout<<summary.str()<<std::flush;
        }
        }
        const auto acceptable=[&]()
        {
            const Eigen::VectorXd residual=actualResidual().cwiseQuotient(scale);
            return x.allFinite() && residual.allFinite() &&
                   residual.stableNorm()<=physicalTolerance*projectedNorm;
        };
        if(!acceptable() && !directFallback && !physicalOperator)
        {
            // Difficult geometries can stagnate with a local preconditioner and
            // restarted GMRES. Retain the complete operator and original residual
            // tolerance; cache a pivoted LU for subsequent RHS at this geometry.
            SolverStageTimer fallback("Climb dense LU fallback factorization (n="+std::to_string(n)+")",true);
            directFallback=std::make_unique<Eigen::PartialPivLU<Eigen::MatrixXd>>(matrix);
            fallback.start("Climb dense LU fallback substitution");
            x=directFallback->solve(rhs);
            residualCurrent=false;
            if(acceptable())
            {
                std::cout<<" (climb direct fallback after "<<iterations<<" GMRES iterations)"<<std::flush;
            }
        }
        if(!acceptable())
        {
            throw failure(physicalOperator?"FMM climb GMRES failed its reduced residual check (no dense fallback).":"Coupled climb solve failed the full reduced residual check after direct fallback.");
        }
        previousSolution=x; // Publish only after the residual check succeeds.
        const Eigen::VectorXd solution=map*(scale.asDiagonal()*x);
        if(profileEnabled)
        {
        const double totalSeconds=std::chrono::duration<double>(ProfileClock::now()-profileStart).count();
        const double otherSeconds=totalSeconds-operatorSeconds-orthogonalizationSeconds-preconditionerSeconds;
        std::ostringstream profile;
        profile<<std::scientific<<std::setprecision(9)
               <<"\n  Climb solve profile: operator="<<(physicalOperator?"FMM":"dense")
               <<", iterations="<<iterations<<", operator_calls="<<operatorCalls
               <<", arnoldi_calls="<<arnoldiCalls<<", residual_calls="<<operatorCalls-arnoldiCalls
               <<", preconditioner_calls="<<preconditionerCalls
               <<", operator_seconds="<<operatorSeconds
               <<", arnoldi_operator_seconds="<<arnoldiOperatorSeconds
               <<", residual_operator_seconds="<<operatorSeconds-arnoldiOperatorSeconds
               <<", orthogonalization_seconds="<<orthogonalizationSeconds
               <<", preconditioner_seconds="<<preconditionerSeconds
               <<", other_seconds="<<otherSeconds<<", total_seconds="<<totalSeconds;
        std::cout<<profile.str()<<std::flush;
        }
        return solution;
    }
}
