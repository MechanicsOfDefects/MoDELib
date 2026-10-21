#ifndef model_ConstrainedVelocitySolve_H_
#define model_ConstrainedVelocitySolve_H_
#include <PeriodicVelocityConstraints.h>
#include <Eigen/SparseQR>
#include <cfloat>
#include <set>

namespace model
{
    using VelocityMatrix = Eigen::SparseMatrix<double>;

    // Solve the existing lumped balance in the admissible velocity space.
    // Reduction must precede solving; projecting unconstrained speeds would
    // not preserve the force balance at shared periodic crossings.
    inline Eigen::VectorXd solveConstrainedVelocityDiagonal(const Eigen::VectorXd &diagonal,
                                                            const Eigen::VectorXd &force, const VelocityMatrix &Z)
    {
        if (diagonal.size() != Z.rows() || force.size() != Z.rows() || !diagonal.allFinite() || !force.allFinite())
        {
            throw std::runtime_error("Invalid constrained velocity balance.");
        }
        if (Z.rows() == Z.cols())
        {
            // The topology eliminator returns identity when there are no constraints.
            Eigen::VectorXd result = Eigen::VectorXd::Zero(force.size());
            for (Eigen::Index k = 0; k < force.size(); ++k)
            {
                if (std::abs(diagonal(k)) > FLT_EPSILON)
                {
                    result(k) = force(k) / diagonal(k);
                }
            }
            return result;
        }
        std::vector<Eigen::Triplet<double>> triplets;
        for (Eigen::Index k = 0; k < diagonal.size(); ++k)
        {
            if (std::abs(diagonal(k)) > FLT_EPSILON)
            {
                triplets.emplace_back(k, k, diagonal(k));
            }
        }
        if (Z.cols() == 0)
        {
            return Eigen::VectorXd::Zero(Z.rows());
        }
        VelocityMatrix D(Z.rows(), Z.rows());
        D.setFromTriplets(triplets.begin(), triplets.end());
        const VelocityMatrix reduced = Z.transpose() * D * Z;
        const Eigen::VectorXd rhs = Z.transpose() * force;
        if (rhs.isZero(0.0))
        {
            return Eigen::VectorXd::Zero(Z.rows());
        }
        Eigen::SparseQR<VelocityMatrix, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(reduced);
        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Constrained climb factorization failed.");
        }
        const Eigen::VectorXd solution = solver.solve(rhs);
        if (!solution.allFinite() || (reduced * solution - rhs).norm() > 1.e-9 * rhs.norm())
        {
            throw std::runtime_error("Inconsistent constrained climb balance.");
        }
        return Z * solution;
    }

    template <int dim, class Network> VelocityMatrix periodicClimbVelocityMap(const Network &network)
    {
        std::vector<Eigen::Matrix<double, dim, 1>> directions;
        for (const auto &entry : network.networkNodes())
        {
            directions.push_back(entry.second.lock()->climbDirection());
        }
        PeriodicVelocityConstraints constraints(directions.size());
        using Vector = Eigen::Matrix<double, dim, 1>;
        using VectorRow = std::map<size_t, Vector>;
        std::map<size_t, VectorRow> interpolants;
        forEachPeriodicVelocityRelation(network,
                                        [&](size_t node, size_t a, size_t b, double wa, double wb)
                                        {
                                            VectorRow expression;
                                            expression.emplace(a, wa * directions[a]);
                                            const auto found = expression.find(b);
                                            if (found == expression.end())
                                            {
                                                expression.emplace(b, wb * directions[b]);
                                            }
                                            else
                                            {
                                                found->second += wb * directions[b];
                                            }
                                            const auto first = interpolants.emplace(node, expression);
                                            if (first.second)
                                            {
                                                // Subdivision tangential motion is a reparameterization, not
                                                // an independent climb DOF. Its source strength uses only the
                                                // projection of the interpolated velocity onto the local normal.
                                                std::map<size_t, double> row{{node, 1.0}};
                                                for (const auto &term : expression)
                                                {
                                                    row[term.first] -= directions[node].dot(term.second);
                                                }
                                                constraints.constrain(std::move(row), false, node);
                                            }
                                            else
                                            {
                                                // Every loop sharing this network node must predict the SAME
                                                // vector motion, not merely the same local normal component.
                                                for (int d = 0; d < dim; ++d)
                                                {
                                                    std::map<size_t, double> row;
                                                    for (const auto &term : first.first->second)
                                                    {
                                                        row[term.first] += term.second(d);
                                                    }
                                                    for (const auto &term : expression)
                                                    {
                                                        row[term.first] -= term.second(d);
                                                    }
                                                    constraints.constrain(std::move(row));
                                                }
                                            }
                                        });
        std::vector<Eigen::Triplet<double>> triplets;
        const auto freeCount = constraints.template append<1>(triplets);
        VelocityMatrix Z(directions.size(), freeCount);
        Z.setFromTriplets(triplets.begin(), triplets.end());
        return Z;
    }

    template <int dim, class Network>
    Eigen::VectorXd interpolateClimbSubdivisionVelocities(const Network &network, const Eigen::VectorXd &nodal)
    {
        Eigen::VectorXd result = nodal;
        std::set<size_t> visited;
        const double tolerance = 1.e-9 * std::max(nodal.norm(), 1.e-30);
        forEachPeriodicVelocityRelation(
            network,
            [&](size_t node, size_t a, size_t b, double wa, double wb)
            {
                const Eigen::Matrix<double, dim, 1> value =
                    wa * nodal.template segment<dim>(dim * a) + wb * nodal.template segment<dim>(dim * b);
                if (visited.insert(node).second)
                {
                    result.template segment<dim>(dim * node) = value;
                }
                else if ((result.template segment<dim>(dim * node) - value).norm() > tolerance)
                {
                    throw std::runtime_error("Incompatible shared climb subdivision velocity.");
                }
            });
        return result;
    }

    template <int dim, class Network>
    Eigen::VectorXd compatiblePredictedVelocities(const Network &network, const Eigen::VectorXd &prediction)
    {
        if (prediction.size() != dim * Eigen::Index(network.networkNodes().size()) || !prediction.allFinite())
        {
            throw std::runtime_error("Invalid Python glide velocity shape or values.");
        }
        std::vector<Eigen::Triplet<double>> triplets;
        const auto freeCount = periodicVelocityConstraints<dim>(network, triplets);
        if (freeCount == network.networkNodes().size())
        {
            return prediction;
        }
        VelocityMatrix Z(prediction.size(), dim * freeCount);
        Z.setFromTriplets(triplets.begin(), triplets.end());
        Eigen::SparseQR<VelocityMatrix, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(Z);
        if (solver.info() != Eigen::Success)
        {
            throw std::runtime_error("Periodic Python velocity projection failed.");
        }
        const Eigen::VectorXd free = solver.solve(prediction);
        if (!free.allFinite())
        {
            throw std::runtime_error("Nonfinite periodic Python velocity projection.");
        }
        return Z * free;
    }
} // namespace model
#endif
