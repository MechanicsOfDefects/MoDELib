/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeIO_H_
#define model_DislocationNodeIO_H_

#include <tuple>
#include <sstream>
#include <ostream>
#include <Eigen/Dense>
#include <iomanip>
#include <Eigen/Dense>

#include <ClusterDynamics/ClusterDynamicsParameters.h>

namespace model
{
    template<int dim> class DislocationNode;
    
    template<short unsigned int dim>
    struct DislocationNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize> ClimbVelocityScalarType;
        
        size_t sID;          // sID
        VectorDim P;          // position
        VectorDim V;          // velocity
        ClimbVelocityScalarType climbVelocityScalar;
        double velocityReduction;             // velocity reduction factor
        int  meshLocation;    // mesh location
//        size_t masterID;          // sID

        
        /**********************************************************************/
        DislocationNodeIO(const DislocationNode<dim>& dn);
        
        /**********************************************************************/
        DislocationNodeIO(const size_t& sID_in,          // sID
                          const VectorDim& P_in,          // position
                          const VectorDim& V_in,          // velocity
                          const ClimbVelocityScalarType& climbVelocityScalar_in,
                          const double& velocityReduction_in,             // velocity reduction factor
//                          const size_t& snID_in,          // component ID
                          const int&  meshLocation_in);
        
        /**********************************************************************/
        DislocationNodeIO();

        /**********************************************************************/
        DislocationNodeIO(std::stringstream& ss);

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationNodeIO<dim>& ds)
        {
            ds.write(os);
            return os;
        }

        std::ostream& write(std::ostream& os) const;
        
	};
	
}
#endif

