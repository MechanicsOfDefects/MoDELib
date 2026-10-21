/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_EllipsoidalInclusion_H_
#define model_EllipsoidalInclusion_H_

#include <memory>
#include <Eigen/Dense>
#include <InclusionMicrostructure/EshelbyInclusionBase.h>

namespace model
{
    template <int dim> class EllipsoidalInclusion : public EshelbyInclusionBase<dim>
    {
        static_assert(dim == 3, "EllipsoidalInclusion requires dim=3");
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        typedef Eigen::Matrix<double, dim, 1> VectorDim;
        typedef Eigen::Matrix<double, 12, 1> Integrals;

      public:
        // Immutable geometry: columns of axes are mutually orthogonal semiaxis vectors.
        const VectorDim C;
        const MatrixDim axes;

        EllipsoidalInclusion(const VectorDim &_C, const MatrixDim &_axes, const MatrixDim &_eT, const double &_nu,
                             const double &_mu, const double &_mobilityReduction, const int &_phaseID,
                             const std::shared_ptr<SecondPhase<dim>> &sph, const double quadratureTolerance = 1.0e-3,
                             const unsigned int quadratureMaxDepth = 30);

        // As in SphericalInclusion, contains excludes the boundary;
        // stress uses the interior value on an exactly classified boundary.
        bool contains(const VectorDim &x) const override;
        MatrixDim stress(const VectorDim &x) const override;

      private:
        double lengthScale;
        MatrixDim rotation;
        VectorDim axisSquared;
        double axisProduct;
        const double tolerance;
        const unsigned int maxDepth;
        Integrals interiorIntegrals;

        struct InteriorStressCache
        {
            MatrixDim eigenstrain;
            MatrixDim stress;
        };
        // Immutable snapshots allow concurrent stress() calls when eT is not modified.
        mutable std::shared_ptr<const InteriorStressCache> interiorStressCache;

        MatrixDim evaluateStress(const VectorDim &r, const bool inside, const MatrixDim &eigenstrain) const;

        VectorDim localPosition(const VectorDim &x) const;
        double ellipsoidalCoordinate(const VectorDim &r) const;
        Integrals integrate(const double ell) const;
    };
} // namespace model
#endif
