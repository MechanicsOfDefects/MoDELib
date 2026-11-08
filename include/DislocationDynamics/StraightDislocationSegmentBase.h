#ifndef model_StraightDislocationSegmentBase_H_
#define model_StraightDislocationSegmentBase_H_
#include <DislocationDynamics/StraightDislocationSegmentMoments.h>
#include <memory>
namespace model
{
    enum class SegmentKernelPart { full, shortRange, longRange };
    struct SegmentCoreParameters
    {
        double coreWidth=1,ewaldWidth=4,mu=1,nu=.3,maxDistance=1;
        int order=8,resolution=64;
    };
    class SegmentConcentrationKernelBase
    {
    public:
        virtual ~SegmentConcentrationKernelBase()=default;
        virtual Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double begin=0,double end=1) const=0;
        virtual Eigen::Vector3d flux(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double,double,double begin=0,double end=1) const=0;
    };
    class StraightDislocationSegmentBase
    {
    public:
        using Segment=helmholtzSegment::Segment;
        using ConcentrationPtr=std::shared_ptr<const SegmentConcentrationKernelBase>;
        virtual ~StraightDislocationSegmentBase()=default;
        virtual const char* name() const=0;
        virtual Eigen::Matrix3d stress(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Matrix3d nearStress(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Matrix3d shortStress(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Matrix3d longStress(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Vector3d displacement(const Eigen::Vector3d&,const Segment&) const=0;
        virtual double elasticInteractionEnergy(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Segment&) const=0;
        virtual ConcentrationPtr concentration(const Eigen::Matrix3d&,bool flux=false,SegmentKernelPart part=SegmentKernelPart::full) const=0;
        virtual SegmentMomentIntegrator::StressMoments stressMoments(const SegmentMomentIntegrator&,const Segment&,const Eigen::Vector3d&,const Eigen::Vector3d&) const=0;
        virtual SegmentMomentIntegrator::ConcentrationMoments concentrationMoments(const SegmentMomentIntegrator&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double,double,const std::array<double,2>&) const=0;
    };
}
#endif
