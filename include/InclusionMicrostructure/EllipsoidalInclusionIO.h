/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */
#ifndef model_EllipsoidalInclusionIO_H_
#define model_EllipsoidalInclusionIO_H_

#include <cstddef>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

namespace model
{
    template <int dim> class EllipsoidalInclusion;
    template <int dim> struct EllipsoidalInclusionIO
    {
        static_assert(dim == 3, "EllipsoidalInclusionIO requires dim=3");
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        typedef Eigen::Matrix<double, dim, 1> VectorDim;

        std::size_t inclusionID;
        VectorDim C;
        MatrixDim axes; // Columns are mutually orthogonal semiaxis vectors.
        MatrixDim eT;
        double mobilityReduction;
        int phaseID;

        EllipsoidalInclusionIO(const EllipsoidalInclusion<dim> &ei);

        EllipsoidalInclusionIO(const std::size_t &inclusionID_in, const VectorDim &C_in, const MatrixDim &axes_in,
                               const MatrixDim &eT_in, const double &mobilityReduction_in, const int &phaseID_in);
        EllipsoidalInclusionIO();
        EllipsoidalInclusionIO(std::stringstream &ss);

        // One record: ID, C, axes (row-major), eT (row-major), mobility, phase.
        // Material constants and numerical settings are supplied when reconstructing
        // the inclusion, as material constants are for SphericalInclusionIO.
        void write(std::ostream &os) const;
        template <class T> friend T &operator<<(T &os, const EllipsoidalInclusionIO<dim> &ei)
        {
            ei.write(os);
            return os;
        }
    };
} // namespace model
#endif
