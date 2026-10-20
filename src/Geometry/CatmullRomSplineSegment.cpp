/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <CatmullRomSplineSegment.h>
#include <math.h>

namespace model
{

    template <int dim>
    auto CatmullRomSplineSegment<dim>::get_M(const VectorDim &P0, const VectorDim &P1, const VectorDim &P2,
                                             const VectorDim &P3) -> Eigen::Matrix<double, 4, 4>
    {
        Eigen::Matrix<double, 4, 4> temp;
        const double g0(std::pow((P1 - P0).norm(), 0.5));
        const double g1(std::pow((P2 - P1).norm(), 0.5));
        const double g2(std::pow((P3 - P2).norm(), 0.5));

        temp << 0, 1.0, 0, 0, -std::pow(g1, 2) / g0 / (g0 + g1), g1 / g0 - 1.0, g0 / (g1 + g0), 0,
            2.0 * std::pow(g1, 2) / g0 / (g0 + g1), -2.0 * g1 / g0 - g1 / (g1 + g2), 2.0 * g1 / (g1 + g0) + g1 / g2,
            -std::pow(g1, 2) / g2 / (g1 + g2), -std::pow(g1, 2) / g0 / (g0 + g1), g1 / g0 + g1 / (g1 + g2),
            -g1 / (g0 + g1) - g1 / g2, std::pow(g1, 2) / g2 / (g1 + g2);

        return temp.transpose();
    }

    template <int dim>
    CatmullRomSplineSegment<dim>::CatmullRomSplineSegment(const VectorDim &P0, const VectorDim &P1, const VectorDim &P2,
                                                          const VectorDim &P3)
        : P((PointMatrixType() << P0, P1, P2, P3).finished()), M(get_M(P0, P1, P2, P3)), PM(P * M)
    {
        // std::cout<<P<<std::endl;
    }

    template <int dim> auto CatmullRomSplineSegment<dim>::get_r(const double &u) const -> VectorDim
    {
        return PM * (Eigen::Matrix<double, 4, 1>() << 1.0, u, std::pow(u, 2), std::pow(u, 3)).finished();
    }

    template <int dim> auto CatmullRomSplineSegment<dim>::get_ru(const double &u) const -> VectorDim
    {
        return PM * (Eigen::Matrix<double, 4, 1>() << 0.0, 1.0, 2.0 * u, 3.0 * std::pow(u, 2)).finished();
    }

    template <int dim> auto CatmullRomSplineSegment<dim>::get_ruu(const double &u) const -> VectorDim
    {
        return PM * (Eigen::Matrix<double, 4, 1>() << 0.0, 0.0, 2.0, 6.0 * u).finished();
    }
    template class CatmullRomSplineSegment<3>;

} // namespace model
