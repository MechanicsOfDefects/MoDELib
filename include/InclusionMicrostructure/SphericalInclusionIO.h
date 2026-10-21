/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionIO_H_
#define model_SphericalInclusionIO_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
// #include <Material.h>


// #include <DislocationStress.h>
//  http://solidmechanics.org/text/Chapter5_4/Chapter5_4.htm


#include <Eigen/Dense>
#include <sstream>
#include <iomanip>

namespace model
{
    template <int dim> class SphericalInclusion;

    template <int dim> struct SphericalInclusionIO
    {
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        typedef Eigen::Matrix<double, dim, 1> VectorDim;

        size_t inclusionID;
        VectorDim C;
        double a; // size (radius)
        MatrixDim eT;
        double mobilityReduction;
        int phaseID;

        /**********************************************************************/
        SphericalInclusionIO(const SphericalInclusion<dim> &ei);

        SphericalInclusionIO(const size_t &inclusionID_in, const VectorDim &C_in, const double &a_in,
                             const MatrixDim &eT_in, const double &mobilityReduction_in, const int &phaseID_in);

        SphericalInclusionIO();

        SphericalInclusionIO(std::stringstream &ss);

        void write(std::ostream &os) const;

        template <class T> friend T &operator<<(T &os, const SphericalInclusionIO<dim> &ei)
        {
            ei.write(os);
            return os;
        }
    };

} // namespace model
#endif
