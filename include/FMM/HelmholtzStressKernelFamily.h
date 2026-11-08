/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_HelmholtzStressKernelFamily_H_
#define model_HelmholtzStressKernelFamily_H_
#include <FMM/FMMKernelFamilyBase.h>
#include <DislocationDynamics/LazarStressParticleKernel.h>
#include <FMM/KokkosStressFMM.h>
namespace model
{
    template<LazarPart Part> struct FMMRegisteredLazarKernel : FMMStressKernelAdapter<LazarStressParticleKernel<Part>>
    {
        using Base=FMMStressKernelAdapter<LazarStressParticleKernel<Part>>;
        FMMRegisteredLazarKernel(double mu,double nu,double a,double h):Base(LazarStressParticleKernel<Part>(mu,nu,a,h)) {}
        KOKKOS_INLINE_FUNCTION double distance(const double* r) const
        { return Kokkos::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); }
    };
    class HelmholtzStressKernelFamily : public FMMKernelFamilyBase<HelmholtzStressKernelFamily>
    {
        FMMRegisteredLazarKernel<LazarPart::full> full_;
        FMMRegisteredLazarKernel<LazarPart::shortRange> short_;
        FMMRegisteredLazarKernel<LazarPart::longRange> long_;
    public:
        static constexpr int sourceComponents=9,fieldComponents=6;
        HelmholtzStressKernelFamily(double mu,double nu,double a,double h,bool split=true):FMMKernelFamilyBase<HelmholtzStressKernelFamily>(split),full_(mu,nu,a,h),short_(mu,nu,a,h),long_(mu,nu,a,h) {}
        auto fullKernel() const { return full_; }
        auto shortKernel() const { return short_; }
        auto longKernel() const { return long_; }
        double searchRadius(double cutoff) const { return cutoff; }
    };
}
#endif
