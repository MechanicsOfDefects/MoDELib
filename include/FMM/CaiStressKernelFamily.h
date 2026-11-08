/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_CaiStressKernelFamily_H_
#define model_CaiStressKernelFamily_H_
#include <DislocationDynamics/CaiStressParticleKernel.h>
namespace model
{
    struct FMMCaiStressKernel : CaiStressParticleKernel
    {
        using CaiStressParticleKernel::CaiStressParticleKernel;
        KOKKOS_INLINE_FUNCTION void add(const double* r,const double* q,double* out) const
        { addSymmetric(r,q,out); }
        KOKKOS_INLINE_FUNCTION double distance(const double* r) const
        { return Kokkos::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); }
    };
    // Full-kernel FMM only. The generic interaction requires all three accessors;
    // no Cai Ewald decomposition is implied by these interface placeholders.
    class CaiStressKernelFamily
    {
        FMMCaiStressKernel kernel_;
    public:
        static constexpr int sourceComponents=9,fieldComponents=6;
        CaiStressKernelFamily(double mu,double nu,double a,double):kernel_(mu,nu,a) {}
        auto fullKernel() const { return kernel_; }
        auto longKernel() const { return kernel_; }
        auto shortKernel() const { return kernel_; }
        double searchRadius(double cutoff) const { return cutoff; }
    };
}
#endif
