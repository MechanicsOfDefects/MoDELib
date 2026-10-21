#ifndef model_ClusterDynamicsBase_H_
#define model_ClusterDynamicsBase_H_
#include <ClusterDynamics/ClusterDynamicsParameters.h>
namespace model
{
    template<int dim>
    struct ClusterDynamicsBase : ClusterDynamicsParameters<dim>
    {
        const bool useClusterDynamicsFEM;
        explicit ClusterDynamicsBase(const DefectiveCrystalBase<dim>& base);
    };
}
#endif
