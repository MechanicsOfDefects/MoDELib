/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationGlideSolverFactory_h_
#define model_DislocationGlideSolverFactory_h_

#include <DislocationVelocitySolverBase.h>
#include <memory>
#include <Eigen/Dense>

namespace model
{
    template<int dim> class DislocationNetwork;
    
    template <int dim>
    struct DislocationGlideSolverBase : public DislocationVelocitySolverBase<dim>
    {
//        const DislocationNetwork<dim>& DN;
        DislocationGlideSolverBase(const DislocationNetwork<dim>& );
//        virtual Eigen::VectorXd getNodeVelocities() const = 0;
    };

    template <int dim>
    struct DislocationGlideSolverFactory
    {
        static std::shared_ptr<DislocationGlideSolverBase<dim>> getGlideSolver(const DislocationNetwork<dim>& DN, const std::string& solverType);
    };
    
}
#endif
