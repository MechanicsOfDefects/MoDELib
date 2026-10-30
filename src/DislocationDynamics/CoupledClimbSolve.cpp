#include <CoupledClimbSolve.h>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
namespace model
{
    void CoupledClimbSolve::compute(const Eigen::MatrixXd& K,const Eigen::SparseMatrix<double>& local,
                                    const Eigen::SparseMatrix<double>& Z)
    {
        ready=false;
        directFallback.reset();
        physicalSize=K.rows();
        if(K.cols()!=physicalSize || Z.rows()!=physicalSize || local.rows()!=physicalSize
           || local.cols()!=physicalSize || !K.allFinite())
        {
            throw std::runtime_error("Invalid coupled climb matrix dimensions or entries.");
        }
        const Eigen::MatrixXd reduced=Z.transpose()*(K*Z);
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
            Eigen::SparseMatrix<double> P=map.transpose()*local*map;
            P=scale.asDiagonal()*P*scale.asDiagonal();
            localSolver.compute(P);
            if(localSolver.info()!=Eigen::Success)
            {
                throw std::runtime_error("Local climb preconditioner factorization failed.");
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
            return Eigen::VectorXd::Zero(physicalSize);
        }
        const Eigen::Index n=matrix.rows();
        const int restart=std::min<Eigen::Index>(50,n);
        const int maxIterations=std::max<Eigen::Index>(500,2*n);
        const double tolerance=1.e-10;
        const double physicalTolerance=1.e-8;
        // ||S^-1 r|| <= ||r|| / min(S). Use this bound for the Arnoldi
        // estimate, but verify both actual residuals before accepting a solve.
        const double estimatedTarget=std::min(tolerance*rhsNorm,
                                             physicalTolerance*scale.minCoeff()*projectedNorm);
        Eigen::VectorXd x=Eigen::VectorXd::Zero(n);
        int iterations=0;
        const auto failure=[&](const char* reason)
        {
            const Eigen::VectorXd r=rhs-matrix*x;
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
            x=directFallback->solve(rhs);
        }
        // Restarted right-preconditioned GMRES, with reorthogonalization.
        while(!directFallback && iterations<maxIterations)
        {
            Eigen::VectorXd residual=rhs-matrix*x;
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
                W.col(j)=localSolver.solve(V.col(j));
                Eigen::VectorXd w=matrix*W.col(j);
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
                if(std::abs(g(j+1))<=estimatedTarget || breakdown)
                {
                    break;
                }
            }
            x+=W.leftCols(count)*H.topLeftCorner(count,count).triangularView<Eigen::Upper>().solve(g.head(count));
        }
        const auto acceptable=[&]()
        {
            const Eigen::VectorXd residual=(rhs-matrix*x).cwiseQuotient(scale);
            return x.allFinite() && residual.allFinite() &&
                   residual.stableNorm()<=physicalTolerance*projectedNorm;
        };
        if(!acceptable() && !directFallback)
        {
            // Difficult geometries can stagnate with a local preconditioner and
            // restarted GMRES. Retain the complete operator and original residual
            // tolerance; cache a pivoted LU for subsequent RHS at this geometry.
            directFallback=std::make_unique<Eigen::PartialPivLU<Eigen::MatrixXd>>(matrix);
            x=directFallback->solve(rhs);
            if(acceptable())
            {
                std::cout<<" (climb direct fallback after "<<iterations<<" GMRES iterations)"<<std::flush;
            }
        }
        if(!acceptable())
        {
            throw failure("Coupled climb solve failed the full reduced residual check after direct fallback.");
        }
        return map*(scale.asDiagonal()*x);
    }
}
