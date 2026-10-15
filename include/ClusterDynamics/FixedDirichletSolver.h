#ifndef model_FixedDirichletSolver_H_
#define model_FixedDirichletSolver_H_

#include<Eigen/SparseCore>
#include<Eigen/SparseLU>
#include<Eigen/IterativeLinearSolvers>

#include<chrono>
#include<iostream>
#include<iomanip>

#include<TrialBase.h>

namespace model
{
    template<typename DirectSolverType,typename IterativeSolverType>
    class FixedDirichletSolver
    {// solves Ax=b given x=T*x1+g. Only b, and g can change
        
        static_assert(std::is_same<typename DirectSolverType::MatrixType,typename IterativeSolverType::MatrixType>::value, "Direct and Iterative MatrixType must be the same");
        typedef typename DirectSolverType::MatrixType SparseMatrixType;

//    public:
//        typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SparseMatrixType;
        
 //   private:

        SparseMatrixType A1;
        
//        Eigen::SparseLU<SparseMatrixType> directSolver;
//        Eigen::BiCGSTAB<SparseMatrixType> iterativeSolver;
        
        DirectSolverType directSolver;
        IterativeSolverType iterativeSolver;

        const std::map<size_t,double>* dirichletConditions;
        const Eigen::Matrix<double,Eigen::Dynamic,1>* dofVector;
        size_t gSize;
        size_t cSize;
        size_t tSize;
        
        SparseMatrixType A;
        SparseMatrixType T;

        public:
        
        const bool use_directSolver;

        const SparseMatrixType& getA() const
        {
            return A;
        }

        const SparseMatrixType& getT() const
        {
            return T;
        }
        
        FixedDirichletSolver(const bool& use_directSolver_in, const double& tol):
        /* init */ dirichletConditions(nullptr)
        /* init */,dofVector(nullptr)
        /* init */,gSize(0)
        /* init */,cSize(0)
        /* init */,tSize(0)
        /* init */,use_directSolver(use_directSolver_in)
        {
            iterativeSolver.setTolerance(tol);
        }
        
        Eigen::VectorXd solve(const Eigen::VectorXd& b) const
        {/*!@param[in] b the rhs of A*x=b
          * @param[in] y the guess for x
          */
            
            Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
            Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
            
            size_t startRow=0;
            size_t col=0;
            
            for (const auto& cIter : *dirichletConditions)
            {
                const size_t& endRow = cIter.first;
                g(endRow)= cIter.second;
                for (size_t row=startRow;row!=endRow;++row)
                {
                    guess(col)=(*dofVector)(row);
                    col++;
                }
                startRow=endRow+1;
            }
            for (size_t row=startRow;row!=gSize;++row)
            {
                guess(col)=(*dofVector)(row);
                col++;
            }
                    
            Eigen::VectorXd b1(T.transpose()*(b-A*g));
            Eigen::VectorXd x(Eigen::VectorXd::Zero(b1.rows()));
            if(use_directSolver)
            {
                x=directSolver.solve(b1);
                if(directSolver.info()!=Eigen::Success)
                {
                    throw std::runtime_error("Direct FixedDirichletSolver failed.");
                }
            }
            else
            {
                x=iterativeSolver.solveWithGuess(b1,guess);
                if(iterativeSolver.info()!=Eigen::Success)
                {
                    throw std::runtime_error("Iterative FixedDirichletSolver failed.");
                }
            }
            return T*x+g;
        }

        static void printSeconds(const std::chrono::steady_clock::time_point& a,const std::chrono::steady_clock::time_point& b)
        {// prints the time from a to b in seconds, with one decimal, and restores the format of std::cout
            const std::ios_base::fmtflags oldFlags(std::cout.flags());
            const std::streamsize oldPrecision(std::cout.precision());
            std::cout<<" "<<std::fixed<<std::setprecision(1)<<std::chrono::duration<double>(b-a).count()<<" s"<<std::flush;
            std::cout.flags(oldFlags);
            std::cout.precision(oldPrecision);
        }

        template<typename BilinearWeakFormType>
//      void compute(const BilinearWeakFormType& bWF)
        void compute(const BilinearWeakFormType& bWF,const bool& printStages=false)
        {
            dirichletConditions=&TrialBase<typename BilinearWeakFormType::TrialFunctionType>::dirichletConditions();
            dofVector=&TrialBase<typename BilinearWeakFormType::TrialFunctionType>::dofVector();
            gSize=TrialBase<typename BilinearWeakFormType::TrialFunctionType>::gSize();
            cSize=dirichletConditions->size();
            tSize = gSize-cSize;

            const auto t0(std::chrono::steady_clock::now());
            if(printStages) std::cout<<"    assembling"<<std::flush;
            A.resize(gSize,gSize);
            const auto aTriplets(bWF.globalTriplets());
            A.setFromTriplets(aTriplets.begin(),aTriplets.end());
            const auto t1(std::chrono::steady_clock::now());
//          if(printStages){printSeconds(t0,t1); std::cout<<" | boundary conditions"<<std::flush;}
            
            std::vector<Eigen::Triplet<double> > tTriplets;
            tTriplets.reserve(tSize);
            size_t startRow=0;
            size_t col=0;
            for (const auto& cIter : *dirichletConditions)
            {
                const size_t& endRow = cIter.first;
                for (size_t row=startRow;row!=endRow;++row)
                {
                    tTriplets.emplace_back(row,col,1.0);
                    col++;
                }
                startRow=endRow+1;
            }
            for (size_t row=startRow;row!=gSize;++row)
            {
                tTriplets.emplace_back(row,col,1.0);
                col++;
            }
            T.resize(gSize,tSize);
            T.setFromTriplets(tTriplets.begin(),tTriplets.end());
            
            A1=T.transpose()*A*T; // store new A1
            const auto t2(std::chrono::steady_clock::now());
//          if(printStages){printSeconds(t1,t2); std::cout<<" | factorizing"<<std::flush;}
            if(printStages){printSeconds(t0,t2); std::cout<<" | factorizing"<<std::flush;}
            
            if(use_directSolver)
            {
                directSolver.compute(A1);
                if(directSolver.info()!=Eigen::Success)
                {
                    throw std::runtime_error("FixedDirichletSolver failed.");
                }
            }
            else
            {
                iterativeSolver.compute(A1);
                if(iterativeSolver.info()!=Eigen::Success)
                {
                    throw std::runtime_error("FixedDirichletSolver failed.");
                }
            }
            const auto t3(std::chrono::steady_clock::now());
//          if(printStages)
//          {// times of the stages, printed for the decomposition of the first step
//              const std::ios_base::fmtflags oldFlags(std::cout.flags());
//              const std::streamsize oldPrecision(std::cout.precision());
//              std::cout<<std::fixed<<std::setprecision(1)
//                       <<"    assembling "<<std::chrono::duration<double>(t1-t0).count()<<" s"
//                       <<" | boundary conditions "<<std::chrono::duration<double>(t2-t1).count()<<" s"
//                       <<" | factorizing "<<std::chrono::duration<double>(t3-t2).count()<<" s"<<std::endl;
//              std::cout.flags(oldFlags);
//              std::cout.precision(oldPrecision);
//          }
            if(printStages){printSeconds(t2,t3); std::cout<<std::endl;}
        }
    };

}
#endif
