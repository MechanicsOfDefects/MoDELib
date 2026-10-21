/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationVelocitySolverBase_h_
#define model_DislocationVelocitySolverBase_h_

#include <memory>
#include <Eigen/Dense>
#include <DislocationDynamicsModule.h>

namespace model
{
    template<int dim> class DislocationNetwork;
    
    template <int dim>
    struct DislocationVelocitySolverBase
    {
        const DislocationNetwork<dim>& DN;

        DislocationVelocitySolverBase(const DislocationNetwork<dim>& );
        virtual ~DislocationVelocitySolverBase() = default;

        virtual Eigen::VectorXd getNodeVelocities() const = 0;
    };
    
}
#endif
