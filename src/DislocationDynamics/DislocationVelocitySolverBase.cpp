/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationVelocitySolverBase_cpp_
#define model_DislocationVelocitySolverBase_cpp_

//#include <PlanarDislocationNode.h>

#include <DislocationVelocitySolverBase.h>

#include <DislocationNetwork.h>

namespace model
{
    
    template <int dim>
    DislocationVelocitySolverBase<dim>::DislocationVelocitySolverBase(const DislocationNetwork<dim>& DN_in) :
    /* init */ DN(DN_in)
    {
        
    }
    
    template struct DislocationVelocitySolverBase<3>;

}
#endif
