/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CatmullRomSplineSegment_H_
#define model_CatmullRomSplineSegment_H_

#include <Eigen/Dense>

#include <ParametricCurve.h>

namespace model
{
    template <int dim> class CatmullRomSplineSegment : public ParametricCurve<CatmullRomSplineSegment<dim>, dim>
    {
        typedef Eigen::Matrix<double, dim, 1> VectorDim;
        typedef Eigen::Matrix<double, dim, 4> PointMatrixType;

        const PointMatrixType P;
        const Eigen::Matrix<double, 4, 4> M;
        const PointMatrixType PM;

        static Eigen::Matrix<double, 4, 4> get_M(const VectorDim &P0, const VectorDim &P1, const VectorDim &P2,
                                                 const VectorDim &P3);

      public:
        CatmullRomSplineSegment(const VectorDim &P0, const VectorDim &P1, const VectorDim &P2, const VectorDim &P3);

        VectorDim get_r(const double &u) const;

        VectorDim get_ru(const double &u) const;

        VectorDim get_ruu(const double &u) const;
    };

} // namespace model
#endif
