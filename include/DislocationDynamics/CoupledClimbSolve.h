#ifndef model_CoupledClimbSolve_H_
#define model_CoupledClimbSolve_H_
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <memory>
namespace model
{
    // Dense physical operator, sparse local right preconditioner. No interaction
    // is dropped from the balance. Constraint reduction precedes preconditioning.
    class CoupledClimbSolve
    {
        Eigen::MatrixXd matrix;
        Eigen::SparseMatrix<double> map;
        Eigen::SparseMatrix<double> inactiveMap;
        Eigen::VectorXd scale;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> localSolver;
        Eigen::Index physicalSize=0;
        bool ready=false;
        mutable std::unique_ptr<Eigen::PartialPivLU<Eigen::MatrixXd>> directFallback;
    public:
        void compute(const Eigen::MatrixXd&,const Eigen::SparseMatrix<double>& local,
                     const Eigen::SparseMatrix<double>& constraints);
        Eigen::VectorXd solve(const Eigen::VectorXd&) const;
        bool usesDirectFallback() const
        {
            return bool(directFallback);
        }
    };
}
#endif
