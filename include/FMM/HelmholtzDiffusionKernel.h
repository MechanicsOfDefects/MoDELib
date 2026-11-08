/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_HelmholtzDiffusionKernel_H_
#define model_HelmholtzDiffusionKernel_H_
#include <DislocationDynamics/CoreType.h>
#include <Kokkos_Core.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <cmath>
#include <stdexcept>
namespace model
{
    enum class HelmholtzDiffusionPart { full, shortRange, longRange };
    // Positive Green function: -div(D grad(g)) = w_a. Weights/signs belong to physics.
    template<HelmholtzDiffusionPart Part> class HelmholtzDiffusionKernel
    {
        CoreType coreType_;
        double metric_[9],a_,h_,d_,ratio_,factor_,series_[16];
    public:
        static constexpr int translationComponents=1;
        HelmholtzDiffusionKernel(const Eigen::Matrix3d& diffusion,double a,double h,CoreType coreType=CoreType::Helmholtz)
            :coreType_(coreType),a_(a),h_(h),ratio_(a/h)
        {
            if(!diffusion.allFinite() || !diffusion.isApprox(diffusion.transpose(),1e-13) ||
               !std::isfinite(a)||!std::isfinite(h)||a<=0||h<=0||ratio_<=0||ratio_>0.25)
            { throw std::invalid_argument("Helmholtz diffusion requires SPD D and h>=4*a>0"); }
            const Eigen::Matrix3d D=(diffusion+diffusion.transpose())*.5;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(D);
            if(es.info()!=Eigen::Success || es.eigenvalues().minCoeff()<=0)
            { throw std::invalid_argument("Diffusion tensor must be positive definite"); }
            d_=std::exp(es.eigenvalues().array().log().mean());
            const Eigen::Matrix3d metric=d_*D.inverse();
            factor_=1/(4*3.14159265358979323846*d_);
            if(!metric.allFinite() || !std::isfinite(factor_/a_) || factor_/a_<=0)
            { throw std::invalid_argument("Diffusion/core scale outside floating point range"); }
            for(int i=0;i<3;++i) {
      for(int j=0;j<3;++j) { metric_[3*i+j]=metric(i,j); }
    }
    double factorial=1;
            for(int n=1;n<=16;++n)
            {
                factorial*=n;
                series_[n-1]=(n%2?1.0:-1.0)*(std::pow(ratio_,n)-ratio_*ratio_)/((1-ratio_*ratio_)*factorial);
            }
        }
        KOKKOS_INLINE_FUNCTION double distance(const double* r) const
        {
            double squared=0;
            for(int i=0;i<3;++i) {
      for(int j=0;j<3;++j) { squared+=r[i]*metric_[3*i+j]*r[j]; }
    }
    return Kokkos::sqrt(squared>0?squared:0);
        }
        KOKKOS_INLINE_FUNCTION double value(const double* r) const
        {
            const double rho=distance(r),z=rho/a_;
            if constexpr(Part==HelmholtzDiffusionPart::full)
            {
                if(coreType_==CoreType::Cai)
                {
                    return factor_/Kokkos::sqrt(rho*rho+a_*a_);
                }
                return rho==0?factor_/a_:factor_*(-Kokkos::expm1(-z))/rho;
            }
            if constexpr(Part==HelmholtzDiffusionPart::shortRange)
            {
                return rho==0?factor_/(a_*(1+ratio_)):
                    factor_*Kokkos::exp(-rho/h_)*(-Kokkos::expm1(-(1-ratio_)*z))/((1-ratio_*ratio_)*rho);
            }
            if(z<.1)
            {
                double v=series_[15]; for(int n=14;n>=0;--n) { v=v*z+series_[n]; }
                return factor_/a_*v;
            }
            return factor_*(-Kokkos::expm1(-rho/h_)+ratio_*ratio_*Kokkos::expm1(-z))/((1-ratio_*ratio_)*rho);
        }
        KOKKOS_INLINE_FUNCTION void add(const double* r,const double* q,double* out) const { out[0]+=value(r)*q[0]; }
        std::array<double,13> translationKey() const
        {
            std::array<double,13> key{};
            for(int i=0;i<9;++i) { key[i]=metric_[i]; }
            key[9]=a_; key[10]=h_; key[11]=d_; key[12]=static_cast<double>(Part)+4*static_cast<double>(coreType_); return key;
        }
        KOKKOS_INLINE_FUNCTION void prepareTranslation(const double* r,double* data) const { data[0]=value(r); }
        KOKKOS_INLINE_FUNCTION void applyTranslation(const double* data,const double* q,double* out) const { out[0]+=data[0]*q[0]; }
    };
    class HelmholtzDiffusionFamily
    {
        HelmholtzDiffusionKernel<HelmholtzDiffusionPart::full> full_;
        HelmholtzDiffusionKernel<HelmholtzDiffusionPart::shortRange> short_;
        HelmholtzDiffusionKernel<HelmholtzDiffusionPart::longRange> long_;
        double radiusScale_;
    public:
        static constexpr int sourceComponents=1,fieldComponents=1;
        HelmholtzDiffusionFamily(const Eigen::Matrix3d& D,double a,double h,CoreType coreType=CoreType::Helmholtz):full_(D,a,h,coreType),short_(D,a,h,coreType),long_(D,a,h,coreType)
        {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es((D+D.transpose())*.5);
            radiusScale_=std::sqrt(es.eigenvalues().maxCoeff()/std::exp(es.eigenvalues().array().log().mean()));
        }
        auto fullKernel() const { return full_; }
        auto shortKernel() const { return short_; }
        auto longKernel() const { return long_; }
        // Cutoff is in the determinant-one diffusion metric, not physical Euclidean distance.
        double searchRadius(double cutoff) const { return cutoff*radiusScale_; }
    };
}
#endif
