#ifdef MODELIB_TEST_EIGEN_FALLBACK
#undef MODELIB_UMFPACK
#endif
#include <WeakProblemSolvers.h>
#include <Eigen/LU>
#include <iostream>
#include <stdexcept>
#include <type_traits>

int main()
{
    using Matrix=Eigen::SparseMatrix<double>;
    using Solver=model::DefaultWeakProblemDirectSolver<Matrix>;
#ifdef MODELIB_UMFPACK
    static_assert(std::is_same_v<Solver,Eigen::UmfPackLU<Matrix>>);
#else
    static_assert(std::is_same_v<Solver,Eigen::SparseLU<Matrix>>);
#endif
    try
    {
        // Nonsymmetric stiffness plus a Lagrange-multiplier border: not SPD.
        Eigen::Matrix4d dense;
        dense<<4,1,0,1, 2,5,1,1, 0,1,3,1, 1,1,1,0;
        Matrix matrix=dense.sparseView();
        Solver solver;
        solver.compute(matrix);
        if(solver.info()!=Eigen::Success)
        {
            throw std::runtime_error("Selected solver cannot factor the bordered system");
        }
        for(double value : {1.0,3.0})
        {
            const Eigen::Vector4d rhs(value,-2,4,0.5);
            const Eigen::VectorXd solution=solver.solve(rhs);
            if(solver.info()!=Eigen::Success || (solution-dense.fullPivLu().solve(rhs)).norm()>1.e-12)
            {
                throw std::runtime_error("Selected solver disagrees with dense reference");
            }
        }
        std::cout<<"PASS solver selection and repeated RHS on an indefinite system\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
