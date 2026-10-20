/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopNodeIO_H_
#define model_DislocationLoopNodeIO_H_

#include <tuple>
#include <sstream>
#include <ostream>
#include <Eigen/Dense>
#include <iomanip>
#include <limits>
#include <Eigen/Dense>


namespace model
{
    template<int dim> class DislocationLoopNode;
    
    template<short unsigned int dim>
    struct DislocationLoopNodeIO
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;

        size_t loopID;
        size_t sID;          // sID
        VectorDim P;          // position
        size_t networkNodeID;          // sID
        VectorDim periodicShift;
        std::pair<short int,short int> edgeIDs;          // sID

        
        /**********************************************************************/
        DislocationLoopNodeIO(const DislocationLoopNode<dim>& dn);
        
        /**********************************************************************/
        DislocationLoopNodeIO(const size_t& sID_in,          // sID
                                const size_t& loopID_in,
                          const VectorDim& P_in,          // position
                          const size_t& nnID,
                            const VectorDim& shift,
                              const  std::pair<short int,short int>& edgeID_in);
        
        /**********************************************************************/
        DislocationLoopNodeIO();

        /**********************************************************************/
        DislocationLoopNodeIO(std::stringstream& ss);

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopNodeIO<dim>& ds)
        {
            ds.write(os);
            return os;
        }

        std::ostream& write(std::ostream& os) const;
        
	};
	
}
#endif
