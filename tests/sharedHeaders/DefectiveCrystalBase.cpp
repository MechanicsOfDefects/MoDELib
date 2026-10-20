#include <DefectiveCrystal/DefectiveCrystalBase.h>

static_assert(sizeof(model::DefectiveCrystalBase<3>)>0);

#include <type_traits>
template<class T,class=void> struct IsComplete : std::false_type {};
template<class T> struct IsComplete<T,std::void_t<decltype(sizeof(T))>> : std::true_type {};
static_assert(!IsComplete<model::DislocationDynamicsBase<3>>::value);
static_assert(!IsComplete<model::ClusterDynamicsBase<3>>::value);
static_assert(!IsComplete<model::ElasticDeformationBase<3>>::value);
