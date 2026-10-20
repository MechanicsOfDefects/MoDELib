/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusionIO_H_
#define model_PolyhedronInclusionIO_H_

#include <tuple>

#include <Eigen/Dense>
#include <sstream>
#include <iomanip>

namespace model
{
    template <int dim> class PolyhedronInclusion;

    template <short unsigned int dim> struct PolyhedronInclusionIO
    {

        typedef Eigen::Matrix<double, dim, dim> MatrixDim;

        size_t inclusionID;
        MatrixDim eT;
        double mobilityReduction;
        int phaseID;

        PolyhedronInclusionIO(const PolyhedronInclusion<dim> &ei);

        PolyhedronInclusionIO(const size_t &inclusionID_in, const MatrixDim &eT_in, const double &mobilityReduction_in,
                              const int &phaseID_in);

        PolyhedronInclusionIO();

        PolyhedronInclusionIO(std::stringstream &ss);

        void write(std::ostream &os) const;

        template <class T> friend T &operator<<(T &os, const PolyhedronInclusionIO<dim> &ei)
        {
            ei.write(os);
            return os;
        }
    };

} // namespace model
#endif
