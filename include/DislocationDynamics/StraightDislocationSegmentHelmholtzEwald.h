#ifndef model_StraightDislocationSegmentHelmholtzEwald_H_
#define model_StraightDislocationSegmentHelmholtzEwald_H_
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
namespace model
{
    class HelmholtzEwaldLongConcentration final : public SegmentConcentrationKernelBase
    {
        std::shared_ptr<const SegmentConcentrationKernelBase> full_,short_;
    public:
        HelmholtzEwaldLongConcentration(std::shared_ptr<const SegmentConcentrationKernelBase> full,std::shared_ptr<const SegmentConcentrationKernelBase> shortPart)
            :full_(std::move(full)),short_(std::move(shortPart)) {}
        Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        {
            return full_->concentrationIntegrals(a,b,x,begin,end)-short_->concentrationIntegrals(a,b,x,begin,end);
        }
        Eigen::Vector3d flux(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        {
            return full_->flux(a,b,x,q0,q1,begin,end)-short_->flux(a,b,x,q0,q1,begin,end);
        }
    };
    class StraightDislocationSegmentHelmholtzEwald final : public StraightDislocationSegmentHelmholtz
    {
        mutable std::once_flag shortOnce_;
        mutable std::unique_ptr<helmholtzSegment::Evaluator> short_;
    public:
        explicit StraightDislocationSegmentHelmholtzEwald(const SegmentCoreParameters& parameters):StraightDislocationSegmentHelmholtz(parameters) {}
        const char* name() const override { return "HelmholtzEwald"; }
        Eigen::Matrix3d nearStress(const Eigen::Vector3d& x,const Segment& s) const override { return shortStress(x,s); }
        Eigen::Matrix3d shortStress(const Eigen::Vector3d& x,const Segment& s) const override
        {
            std::call_once(shortOnce_,[this]()
            {
                const auto& p=parameters_;
                short_=std::make_unique<helmholtzSegment::Evaluator>(p.coreWidth,p.ewaldWidth,p.mu,p.nu,p.order,p.resolution,p.maxDistance,false);
            });
            return short_->stress(x,s);
        }
        Eigen::Matrix3d longStress(const Eigen::Vector3d& x,const Segment& s) const override { return stress(x,s)-shortStress(x,s); }
        ConcentrationPtr concentration(const Eigen::Matrix3d& D,bool flux=false,SegmentKernelPart part=SegmentKernelPart::full) const override
        {
            if(part==SegmentKernelPart::full)
            {
                return StraightDislocationSegmentHelmholtz::concentration(D,flux);
            }
            auto shortPart=std::make_shared<HelmholtzSegmentConcentration>(D,parameters_,true,flux);
            if(part==SegmentKernelPart::shortRange)
            {
                return shortPart;
            }
            return std::make_shared<HelmholtzEwaldLongConcentration>(StraightDislocationSegmentHelmholtz::concentration(D,flux),shortPart);
        }
    };
}
#endif
