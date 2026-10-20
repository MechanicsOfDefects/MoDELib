#include <ClusterDynamics/ClusterDynamicsBase.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
namespace model
{
    template<int dim>
    ClusterDynamicsBase<dim>::ClusterDynamicsBase(const DefectiveCrystalBase<dim>& base) :
        ClusterDynamicsParameters<dim>(base),
        useClusterDynamicsFEM(base.fe && TextFileParser(base.simulationParameters.traitsIO.inputFilesFolder+"/ClusterDynamics.txt").readScalar<int>("useClusterDynamicsFEM",true))
    {
    }
    template struct ClusterDynamicsBase<3>;
}
