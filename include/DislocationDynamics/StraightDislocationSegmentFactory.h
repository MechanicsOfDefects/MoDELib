#ifndef model_StraightDislocationSegmentFactory_H_
#define model_StraightDislocationSegmentFactory_H_
#include <DislocationDynamics/CoreType.h>
#include <DislocationDynamics/StraightDislocationSegmentCai.h>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtzEwald.h>
namespace model
{
    inline std::shared_ptr<const StraightDislocationSegmentBase> makeStraightDislocationSegmentCore(CoreType type,bool ewald,const SegmentCoreParameters& parameters)
    {
        switch(type)
        {
            case CoreType::Cai:
                if(ewald)
                {
                    throw std::invalid_argument("Cai Ewald screening is not implemented");
                }
                return std::make_shared<StraightDislocationSegmentCai>(parameters);
            case CoreType::Helmholtz:
                if(ewald)
                {
                    return std::make_shared<StraightDislocationSegmentHelmholtzEwald>(parameters);
                }
                return std::make_shared<StraightDislocationSegmentHelmholtz>(parameters);
        }
        throw std::invalid_argument("Unknown segment core type");
    }
}
#endif
