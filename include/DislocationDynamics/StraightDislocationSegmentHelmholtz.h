#ifndef model_StraightDislocationSegmentHelmholtz_H_
#define model_StraightDislocationSegmentHelmholtz_H_
#include <DislocationDynamics/StraightDislocationSegmentBase.h>
#include <DislocationDynamics/HelmholtzSegmentTables.h>
namespace model
{
    class HelmholtzSegmentConcentration final : public SegmentConcentrationKernelBase
    {
        helmholtzSegment::ConcentrationEvaluator table_;
    public:
        HelmholtzSegmentConcentration(const Eigen::Matrix3d& D,const SegmentCoreParameters& p,bool screened,bool flux)
            :table_(D,p.coreWidth,p.ewaldWidth,screened,p.order,p.resolution,p.maxDistance,flux) {}
        Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        {
            return table_.concentrationIntegrals(a,b,x,begin,end);
        }
        Eigen::Vector3d flux(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        {
            return table_.flux(a,b,x,q0,q1,begin,end);
        }
    };
    class StraightDislocationSegmentHelmholtz : public StraightDislocationSegmentBase
    {
    protected:
        SegmentCoreParameters parameters_;
    private:
        mutable std::once_flag stressOnce_;
        mutable std::unique_ptr<helmholtzSegment::Evaluator> full_;
    public:
        explicit StraightDislocationSegmentHelmholtz(const SegmentCoreParameters& parameters):parameters_(parameters) {}
        const char* name() const override { return "Helmholtz"; }
        Eigen::Matrix3d stress(const Eigen::Vector3d& x,const Segment& s) const override
        {
            std::call_once(stressOnce_,[this]()
            {
                const auto& p=parameters_;
                full_=std::make_unique<helmholtzSegment::Evaluator>(p.coreWidth,p.ewaldWidth,p.mu,p.nu,p.order,p.resolution,p.maxDistance,true);
            });
            return full_->stress(x,s);
        }
        Eigen::Matrix3d nearStress(const Eigen::Vector3d& x,const Segment& s) const override { return stress(x,s); }
        Eigen::Matrix3d shortStress(const Eigen::Vector3d&,const Segment&) const override { throw std::logic_error("Select HelmholtzEwald for screened kernels"); }
        Eigen::Matrix3d longStress(const Eigen::Vector3d&,const Segment&) const override { throw std::logic_error("Select HelmholtzEwald for screened kernels"); }
        Eigen::Vector3d displacement(const Eigen::Vector3d&,const Segment&) const override { throw std::logic_error("Helmholtz displacement is not implemented"); }
        double elasticInteractionEnergy(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Segment&) const override { throw std::logic_error("Helmholtz interaction energy is not implemented"); }
        ConcentrationPtr concentration(const Eigen::Matrix3d& D,bool flux=false,SegmentKernelPart part=SegmentKernelPart::full) const override
        {
            if(part!=SegmentKernelPart::full)
            {
                throw std::logic_error("Select HelmholtzEwald for screened kernels");
            }
            return std::make_shared<HelmholtzSegmentConcentration>(D,parameters_,false,flux);
        }
        SegmentMomentIntegrator::StressMoments stressMoments(const SegmentMomentIntegrator& helper,const Segment& segment,const Eigen::Vector3d& center,const Eigen::Vector3d& width) const override
        {
            return helper.stress(segment,center,width);
        }
        SegmentMomentIntegrator::ConcentrationMoments concentrationMoments(const SegmentMomentIntegrator& helper,const Eigen::Vector3d& p0,const Eigen::Vector3d& p1,const Eigen::Vector3d& center,const Eigen::Vector3d& width,double begin,double end,const std::array<double,2>& weights) const override
        {
            return helper.concentration(p0,p1,center,width,begin,end,weights);
        }
    };
}
#endif
