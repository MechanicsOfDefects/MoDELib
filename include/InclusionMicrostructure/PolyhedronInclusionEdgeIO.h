/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusionEdgeIO_H_
#define model_PolyhedronInclusionEdgeIO_H_

#include <tuple>

#include <Eigen/Dense>
#include <sstream>
#include <iomanip>

namespace model
{

    struct PolyhedronInclusionEdgeIO
    {

        //        typedef Eigen::Matrix<double,dim,1> VectorDim;

        size_t inclusionID; // sID
        size_t faceID;      // sID
        size_t sourceID;    // sID
        size_t sinkID;      // sID

        /**********************************************************************/
        PolyhedronInclusionEdgeIO(const size_t &inclusionID_in, const size_t &faceID_in, const size_t &sourceID_in,
                                  const size_t &sinkID_in);

        /**********************************************************************/
        PolyhedronInclusionEdgeIO();

        /**********************************************************************/
        PolyhedronInclusionEdgeIO(std::stringstream &ss);

        /**********************************************************************/
        void write(std::ostream &os) const;

        template <class T> friend T &operator<<(T &os, const PolyhedronInclusionEdgeIO &ds)
        {
            ds.write(os);
            return os;
        }
    };

} // namespace model
#endif
