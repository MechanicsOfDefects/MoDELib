#ifndef model_StraightDislocationSegmentHelmholtzEwald_H_
#define model_StraightDislocationSegmentHelmholtzEwald_H_
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
namespace model
{
    class HelmholtzEwaldConcentration final : public SegmentConcentrationKernelBase
    {
        Eigen::Matrix3d diffusion_;
        SegmentCoreParameters parameters_;
        bool includeFlux_;
        mutable std::once_flag fullOnce_,nearOnce_;
        mutable std::unique_ptr<HelmholtzSegmentConcentration> full_,near_;
        const HelmholtzSegmentConcentration& fullTable() const
        {
            std::call_once(fullOnce_,[this]() { full_=std::make_unique<HelmholtzSegmentConcentration>(diffusion_,parameters_,false,includeFlux_); });
            return *full_;
        }
        const HelmholtzSegmentConcentration& nearTable() const
        {
            std::call_once(nearOnce_,[this]() { near_=std::make_unique<HelmholtzSegmentConcentration>(diffusion_,parameters_,true,includeFlux_); });
            return *near_;
        }
    public:
        HelmholtzEwaldConcentration(const Eigen::Matrix3d& D,const SegmentCoreParameters& p,bool flux)
            :SegmentConcentrationKernelBase(true),diffusion_(D),parameters_(p),includeFlux_(flux) {}
        Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        { return fullTable().concentrationIntegrals(a,b,x,begin,end); }
        Eigen::Vector2d concentrationIntegralsNear(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        { return nearTable().concentrationIntegrals(a,b,x,begin,end); }
        Eigen::Vector2d concentrationIntegralsFar(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        { return concentrationIntegrals(a,b,x,begin,end)-concentrationIntegralsNear(a,b,x,begin,end); }
        Eigen::Vector3d flux(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        { return fullTable().flux(a,b,x,q0,q1,begin,end); }
        Eigen::Vector3d fluxNear(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        { return nearTable().flux(a,b,x,q0,q1,begin,end); }
        Eigen::Vector3d fluxFar(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        { return flux(a,b,x,q0,q1,begin,end)-fluxNear(a,b,x,q0,q1,begin,end); }
    };
    class StraightDislocationSegmentHelmholtzEwald final : public StraightDislocationSegmentHelmholtz
    {
        mutable std::once_flag shortOnce_;
        mutable std::unique_ptr<helmholtzSegment::Evaluator> short_;
    public:
        explicit StraightDislocationSegmentHelmholtzEwald(const SegmentCoreParameters& parameters):StraightDislocationSegmentHelmholtz(parameters,true) {}
        const char* name() const override { return "HelmholtzEwald"; }
        Eigen::Matrix3d stressNear(const Eigen::Vector3d& x,const Segment& s) const override
        {
            std::call_once(shortOnce_,[this]()
            {
                const auto& p=parameters_;
                short_=std::make_unique<helmholtzSegment::Evaluator>(p.coreWidth,p.ewaldWidth,p.mu,p.nu,p.order,p.resolution,p.maxDistance,false);
            });
            return short_->stress(x,s);
        }
        Eigen::Matrix3d stressFar(const Eigen::Vector3d& x,const Segment& s) const override { return stress(x,s)-stressNear(x,s); }
        ConcentrationPtr concentration(const Eigen::Matrix3d& D,bool flux=false) const override
        {
            return std::make_shared<HelmholtzEwaldConcentration>(D,parameters_,flux);
        }
    };
}
#endif
