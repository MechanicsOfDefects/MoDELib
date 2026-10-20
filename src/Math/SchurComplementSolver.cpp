/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <SchurComplementSolver.h>
#include <assert.h>

namespace model
{

    auto SchurComplementSolver::checkDecomposition() const -> void
    {
        switch (llt.info())
        {
        case Eigen::Success:
            break;

        case Eigen::NumericalIssue:
            assert(0 && "SchurComplementSolver found Eigen::NumericalIssue.");
            break;

        case Eigen::NoConvergence:
            assert(0 && "SchurComplementSolver found Eigen::NoConvergence.");
            break;

        case Eigen::InvalidInput:
            assert(0 && "SchurComplementSolver found Eigen::InvalidInput.");
            break;

        default:
            assert(0 && "SchurComplementSolver failed.");
            break;
        }
    }

    SchurComplementSolver::SchurComplementSolver(const Eigen::MatrixXd &KQQ, const Eigen::MatrixXd &KPQ,
                                                 const Eigen::VectorXd &Fq, const Eigen::VectorXd &Fp)
        : llt(KQQ), pplu(KPQ * llt.solve(KPQ.transpose()))
    {
        checkDecomposition();
        // solve(KPQ,Fq,Fp);
        L = pplu.solve(KPQ * llt.solve(Fq) - Fp);
        X = llt.solve(Fq - KPQ.transpose() * L);
    }

    SchurComplementSolver::SchurComplementSolver(const Eigen::MatrixXd &KQQ, const Eigen::MatrixXd &KPQ,
                                                 const Eigen::VectorXd &Fq)
        : llt(KQQ), pplu(KPQ * llt.solve(KPQ.transpose()))
    {
        checkDecomposition();
        // solve(KPQ,Fq);
        L = pplu.solve(KPQ * llt.solve(Fq));
        X = llt.solve(Fq - KPQ.transpose() * L);
    }
} // namespace model
