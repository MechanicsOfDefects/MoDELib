#ifndef model_WeakProblemSolvers_H_
#define model_WeakProblemSolvers_H_

#include <Eigen/Sparse>
#ifdef MODELIB_UMFPACK
#include <Eigen/UmfPackSupport>
#endif

namespace model
{
    // Global constraints and rigid-motion penalties produce indefinite bordered
    // systems. Use general sparse LU, including when CHOLMOD is also available.
    template<class MatrixType>
#ifdef MODELIB_UMFPACK
    using DefaultWeakProblemDirectSolver=Eigen::UmfPackLU<MatrixType>;
#else
    using DefaultWeakProblemDirectSolver=Eigen::SparseLU<MatrixType>;
#endif
}
#endif
