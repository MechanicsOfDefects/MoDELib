/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusionNodeIO_H_
#define model_PolyhedronInclusionNodeIO_H_

#include <tuple>

#include <Eigen/Dense>
#include <sstream>
#include <iomanip>

namespace model
{

    template <short unsigned int dim> struct PolyhedronInclusionNodeIO
    {

        typedef Eigen::Matrix<double, dim, 1> VectorDim;

        //        size_t inclusionID;
        size_t nodeID; // sID
        VectorDim P;

        /**********************************************************************/
        PolyhedronInclusionNodeIO(

            const size_t &nodeID_in, const VectorDim &P_in);

        /**********************************************************************/
        PolyhedronInclusionNodeIO();

        /**********************************************************************/
        PolyhedronInclusionNodeIO(std::stringstream &ss);

        /**********************************************************************/
        void write(std::ostream &os) const;

        template <class T> friend T &operator<<(T &os, const PolyhedronInclusionNodeIO<dim> &ds)
        {
            ds.write(os);
            return os;
        }
    };

} // namespace model
#endif
