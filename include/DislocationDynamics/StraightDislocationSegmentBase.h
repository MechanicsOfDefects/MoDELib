#ifndef model_StraightDislocationSegmentBase_H_
#define model_StraightDislocationSegmentBase_H_
#include <DislocationDynamics/StraightDislocationSegmentMoments.h>
#include <memory>
namespace model
{
    struct SegmentCoreParameters
    {
        double coreWidth=1,ewaldWidth=4,mu=1,nu=.3,maxDistance=1;
        int order=8,resolution=64;
    };
    class SegmentConcentrationKernelBase
    {
    public:
        const bool isSplit;
        explicit SegmentConcentrationKernelBase(bool split=false):isSplit(split) {}
        virtual ~SegmentConcentrationKernelBase()=default;
        virtual Eigen::Vector2d concentrationIntegrals(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double begin=0,double end=1) const=0;
        virtual Eigen::Vector2d concentrationIntegralsNear(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double begin=0,double end=1) const=0;
        virtual Eigen::Vector2d concentrationIntegralsFar(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double begin=0,double end=1) const=0;
        virtual Eigen::Vector3d flux(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double,double,double begin=0,double end=1) const=0;
        virtual Eigen::Vector3d fluxNear(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double,double,double begin=0,double end=1) const=0;
        virtual Eigen::Vector3d fluxFar(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double,double,double begin=0,double end=1) const=0;
    };
    class UnsplitSegmentConcentrationKernel : public SegmentConcentrationKernelBase
    {
    public:
        Eigen::Vector2d concentrationIntegralsNear(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        { return concentrationIntegrals(a,b,x,begin,end); }
        Eigen::Vector2d concentrationIntegralsFar(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double begin=0,double end=1) const override
        { return concentrationIntegrals(a,b,x,begin,end); }
        Eigen::Vector3d fluxNear(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        { return flux(a,b,x,q0,q1,begin,end); }
        Eigen::Vector3d fluxFar(const Eigen::Vector3d& a,const Eigen::Vector3d& b,const Eigen::Vector3d& x,double q0,double q1,double begin=0,double end=1) const override
        { return flux(a,b,x,q0,q1,begin,end); }
    };
    class StraightDislocationSegmentBase
    {
    public:
        using Segment=helmholtzSegment::Segment;
        using ConcentrationPtr=std::shared_ptr<const SegmentConcentrationKernelBase>;
        const bool isSplit;
        explicit StraightDislocationSegmentBase(bool split=false):isSplit(split) {}
        virtual ~StraightDislocationSegmentBase()=default;
        virtual const char* name() const=0;
        virtual Eigen::Matrix3d stress(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Matrix3d stressNear(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Matrix3d stressFar(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Vector3d displacement(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Vector3d displacementNear(const Eigen::Vector3d&,const Segment&) const=0;
        virtual Eigen::Vector3d displacementFar(const Eigen::Vector3d&,const Segment&) const=0;
        virtual double elasticInteractionEnergy(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Segment&) const=0;
        virtual double elasticInteractionEnergyNear(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Segment&) const=0;
        virtual double elasticInteractionEnergyFar(const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Segment&) const=0;
        virtual ConcentrationPtr concentration(const Eigen::Matrix3d&,bool flux=false) const=0;
        virtual SegmentMomentIntegrator::StressMoments stressMoments(const SegmentMomentIntegrator&,const Segment&,const Eigen::Vector3d&,const Eigen::Vector3d&) const=0;
        virtual SegmentMomentIntegrator::ConcentrationMoments concentrationMoments(const SegmentMomentIntegrator&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,const Eigen::Vector3d&,double,double,const std::array<double,2>&) const=0;
    };
}
#endif
