#ifndef model_ElasticDeformationBase_H_
#define model_ElasticDeformationBase_H_
namespace model
{
    template<int dim> struct DefectiveCrystalBase;
    template<int dim>
    struct ElasticDeformationBase
    {
        const bool useElasticDeformationFEM;
        const double inertiaReliefPenaltyFactor;
        explicit ElasticDeformationBase(const DefectiveCrystalBase<dim>& base);
    };
}
#endif
