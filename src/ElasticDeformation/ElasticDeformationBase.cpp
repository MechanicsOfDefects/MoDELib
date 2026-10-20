#include <ElasticDeformation/ElasticDeformationBase.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
namespace model
{
    template<int dim>
    ElasticDeformationBase<dim>::ElasticDeformationBase(const DefectiveCrystalBase<dim>& base) :
        useElasticDeformationFEM(base.fe && TextFileParser(base.simulationParameters.traitsIO.inputFilesFolder+"/ElasticDeformation.txt").readScalar<int>("useElasticDeformationFEM",true)),
        inertiaReliefPenaltyFactor(useElasticDeformationFEM ? TextFileParser(base.simulationParameters.traitsIO.inputFilesFolder+"/ElasticDeformation.txt").readScalar<double>("inertiaReliefPenaltyFactor",true) : 0.0)
    {
    }
    template struct ElasticDeformationBase<3>;
}
