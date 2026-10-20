/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationClimbSolverFactory_h_
#define model_DislocationClimbSolverFactory_h_

#include <DislocationVelocitySolverBase.h>
#include <memory>
#include <Eigen/Dense>
#include <ClusterDynamics.h>
#include <TypeTraits.h>
#include <ClusterDynamics/ClusterDynamicsParameters.h>


namespace model
{
    template<int dim> class DislocationNetwork;
    
    template <int dim>
    struct DislocationClimbSolverBase : public DislocationVelocitySolverBase<dim>
    /*                               */,private std::vector<Eigen::Array<double,1,ClusterDynamicsParameters<TypeTraits<DislocationNetwork<dim>>::dim>::mSize>>

    {
//        const DislocationNetwork<dim>& DN;
        typedef std::vector<Eigen::Array<double,1,ClusterDynamicsParameters<dim>::mSize>> ScalarVelocitiesContainerType;
        const double glideEquilibriumRate;
        const ClusterDynamics<dim>* const CD;
        const double vClimbRef;
        
        DislocationClimbSolverBase(const DislocationNetwork<dim>&,const ClusterDynamics<dim>* const);
//        const ClusterDynamics<dim>*  getCD() const;
        double getVclimbRef() const;
        virtual void computeClimbScalarVelocities(const bool& updateStiffness) =0;

        const ScalarVelocitiesContainerType& scalarVelocities() const;
        ScalarVelocitiesContainerType& scalarVelocities();

        
  //      virtual Eigen::VectorXd getNodeVelocities() const = 0;
    };

    template <int dim>
    struct DislocationClimbSolverFactory
    {

//        static const ClusterDynamics<dim>*  getCD(const DislocationNetwork<dim>& DN);

        static std::shared_ptr<DislocationClimbSolverBase<dim>> getClimbSolver(const DislocationNetwork<dim>& DN, const std::string& solverType);
    };
    
}
#endif
