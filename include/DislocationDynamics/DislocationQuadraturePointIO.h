/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationQuadraturePointIO_H_
#define model_DislocationQuadraturePointIO_H_

#include <Eigen/Dense>
#include <sstream>
#include <ostream>
#include <ClusterDynamics/ClusterDynamicsParameters.h>


namespace model
{
    template<int dim> struct DislocationQuadraturePoint;
    template<int dim>
    struct DislocationQuadraturePointIO
    {

        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        size_t sourceID;
        size_t sinkID;
        int qID;
        VectorDim r;                          // position
        double j;                             // jacobian dl/du
        VectorDim rl;                         // unit tangent dr/dl
        
        MatrixDim stress;
        VectorDim pkForce;
        VectorDim stackingFaultForce;
        VectorDim lineTensionForce;
        VectorDim velocity;
        double elasticEnergyPerLength;
        double coreEnergyPerLength;
        Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize> cCD;
        Eigen::Matrix<double,1,ClusterDynamicsParameters<dim>::mSize> cDD;

        /**********************************************************************/
        DislocationQuadraturePointIO();
        
        /**********************************************************************/
        DislocationQuadraturePointIO(const DislocationQuadraturePoint<dim>& qPoint);
        
        /**********************************************************************/
        DislocationQuadraturePointIO(std::stringstream& ss);
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationQuadraturePointIO<dim>& ds)
        {
            ds.write(os);
            return os;
        }

        std::ostream& write(std::ostream& os) const;
        
    };
    
}
#endif
