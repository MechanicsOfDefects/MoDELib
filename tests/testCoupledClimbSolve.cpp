#include <DislocationDynamics/CoupledClimbSolve.h>
#include <iostream>
#include <stdexcept>
#include <cmath>
using namespace model;
void require(bool ok,const char* message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
int main()
{
    try
    {
        const int n=80;
        Eigen::MatrixXd K(n,n);
        for(int i=0;i<n;++i)
        {
            for(int j=0;j<n;++j)
            {
                K(i,j)=1.0/std::sqrt(0.25+double((i-j)*(i-j)));
                // Source-dependent bias: deliberately nonsymmetric.
                K(i,j)*=1.0+0.003*j;
            }
        }
        Eigen::SparseMatrix<double> P(n,n),Z(n,n-2);
        std::vector<Eigen::Triplet<double>> entries,relations;
        for(int i=0;i<n;++i)
        {
            for(int j=std::max(0,i-1);j<=std::min(n-1,i+1);++j)
            {
                entries.emplace_back(i,j,K(i,j));
            }
            if(i<n-2)
            {
                relations.emplace_back(i,i,1.0);
            }
        }
        // A periodic duplicate and a corner interpolation with opposite direction.
        relations.emplace_back(n-2,0,1.0);
        relations.emplace_back(n-1,1,0.3);
        relations.emplace_back(n-1,2,-0.7);
        P.setFromTriplets(entries.begin(),entries.end());
        Z.setFromTriplets(relations.begin(),relations.end());
        CoupledClimbSolve solver;
        solver.compute(K,P,Z);
        for(int rhs=0;rhs<3;++rhs)
        {
            Eigen::VectorXd f(n);
            for(int i=0;i<n;++i)
            {
                f(i)=std::sin(0.3*i+rhs);
            }
            const Eigen::MatrixXd reduced=Z.transpose()*(K*Z);
            const Eigen::VectorXd exact=Z*reduced.partialPivLu().solve(Z.transpose()*f);
            const Eigen::VectorXd v=solver.solve(f);
            require((v-exact).norm()<1.e-8*exact.norm(),"Dense constrained reference mismatch");
            require((Z.transpose()*(K*v-f)).norm()<1.e-8*(Z.transpose()*f).norm(),"Full residual mismatch");
        }
        require(solver.solve(Eigen::VectorXd::Zero(n)).isZero(0),"Zero drive must give zero velocity");
        // Opposite signed source rates: positive driving terms at both sites
        // can still require absorption at the second site through coupling.
        Eigen::MatrixXd pair(2,2);
        pair<<2,1,1,2;
        Eigen::SparseMatrix<double> identity(2,2);
        identity.setIdentity();
        solver.compute(pair,2.0*identity,identity);
        Eigen::Vector2d force(1,0.1);
        const Eigen::VectorXd speeds=solver.solve(force);
        require(speeds(0)>0 && speeds(1)<0,"Emission must induce opposite-signed absorption");
        require((force.array()/pair.rowwise().sum().array()).minCoeff()>0,"Fixture must expose lumping error");
        // A small scaled residual can conceal an order-one physical residual.
        Eigen::Matrix2d disparate;
        disparate<<1.0,1.0,1.0,1.e24;
        Eigen::SparseMatrix<double> diagonal(2,2);
        diagonal.insert(0,0)=1.0;
        diagonal.insert(1,1)=1.e24;
        solver.compute(disparate,diagonal,identity);
        for(const double amplitude : {1.0,1.e-23,1.e-180})
        {
            const Eigen::Vector2d drive=amplitude*Eigen::Vector2d::Ones();
            const Eigen::VectorXd velocity=solver.solve(drive);
            require((disparate*velocity-drive).stableNorm()<1.e-8*drive.stableNorm(),"Scaled convergence concealed a physical residual");
            require((velocity-Eigen::Vector2d(amplitude,0)).stableNorm()<1.e-8*drive.stableNorm(),"Uneven-scale solve differs from exact solution");
        }
        // Broad spectrum causes restarted GMRES with a diagonal local
        // preconditioner to stagnate; a cached full LU must recover the solution.
        const int hardSize=96;
        Eigen::MatrixXd basis(hardSize,hardSize);
        Eigen::VectorXd eigenvalues(hardSize);
        for(int i=0;i<hardSize;++i)
        {
            eigenvalues(i)=std::pow(10.0,-9.0*i/(hardSize-1));
            for(int j=0;j<hardSize;++j)
            {
                basis(i,j)=std::sqrt(2.0/(hardSize+1))*std::sin(std::acos(-1.0)*(i+1)*(j+1)/(hardSize+1));
            }
        }
        const Eigen::MatrixXd hard=basis*eigenvalues.asDiagonal()*basis.transpose();
        Eigen::SparseMatrix<double> hardIdentity(hardSize,hardSize),hardLocal(hardSize,hardSize);
        hardIdentity.setIdentity();
        for(int i=0;i<hardSize;++i)
        {
            hardLocal.insert(i,i)=hard(i,i);
        }
        solver.compute(hard,hardLocal,hardIdentity);
        require(!solver.usesDirectFallback(),"New geometry retained a stale direct factorization");
        for(int load=0;load<2;++load)
        {
            Eigen::VectorXd drive(hardSize);
            for(int i=0;i<hardSize;++i)
            {
                drive(i)=std::sin(0.31*i+load);
            }
            const Eigen::VectorXd velocity=solver.solve(drive);
            require((hard*velocity-drive).norm()<1.e-8*drive.norm(),"Direct fallback failed physical residual");
            require(solver.usesDirectFallback(),"Stagnating solve did not use cached direct fallback");
        }
        Eigen::MatrixXd inactive=Eigen::MatrixXd::Zero(3,3);
        inactive.topLeftCorner(2,2)=pair;
        Eigen::SparseMatrix<double> I3(3,3);
        I3.setIdentity();
        solver.compute(inactive,inactive.sparseView(),I3);
        require(!solver.usesDirectFallback(),"Matrix rebuild did not invalidate direct fallback");
        Eigen::Vector3d f3(1,0.1,0);
        require(solver.solve(f3)(2)==0,"Inactive node must stay stationary");
        f3(2)=1;
        bool rejected=false;
        try
        {
            solver.solve(f3);
        }
        catch(const std::runtime_error&)
        {
            rejected=true;
        }
        require(rejected,"Inconsistent inactive drive must be rejected");
        std::cout<<"PASS: coupled climb, local preconditioning, mixed signs, constraints, cached RHS and inactive nodes\n";
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
}
