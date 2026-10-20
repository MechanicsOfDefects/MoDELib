/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopLinkIO_h_
#define model_DislocationLoopLinkIO_h_

#include <iomanip>
#include <tuple>
#include <sstream>
#include <ostream>
#include <Eigen/Dense>
#include <TypeTraits.h>

namespace model
{
    template<int dim> class DislocationLoopLink;
    
    template<short unsigned int dim>
    struct DislocationLoopLinkIO
    {
        
        
        size_t loopID;          // sID
        size_t sourceID;          // sID
        size_t sinkID;          // sID
        bool hasNetworkLink;
        int  meshLocation;    // mesh location
        
        /**********************************************************************/
        DislocationLoopLinkIO(const DislocationLoopLink<dim>& ll);
        
        /**********************************************************************/
        DislocationLoopLinkIO(const size_t& loopID_in,          // sID
                          const size_t& sourceID_in,          // position
                          const size_t& sinkID_in,          // velocity
                              const bool& hasNetworkLink_in,
                          const int& meshLocation_in);
        
        /**********************************************************************/
        DislocationLoopLinkIO();
        
        /**********************************************************************/
        DislocationLoopLinkIO(std::stringstream& ss);

        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopLinkIO<dim>& ds)
        {
            ds.write(os);
            return os;
        }

        std::ostream& write(std::ostream& os) const;
        
	};
	
}
#endif

