/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosStressFMM_H_
#define model_KokkosStressFMM_H_
#include <FMM/KokkosParticleFMM.h>
namespace model
{
    // Adapter preserves legacy stress kernels, including their seven-coefficient cache.
    template<class Kernel> struct FMMStressKernelAdapter : Kernel
    {
        static constexpr int translationComponents=7;
        explicit FMMStressKernelAdapter(const Kernel& k):Kernel(k) {}
        KOKKOS_INLINE_FUNCTION void add(const double* r,const double* q,double* out) const
        { this->addSymmetric(r,q,out); }
    };
    struct FMMStressMatrixOutput
    {
        static constexpr int components=9;
        template<class View> KOKKOS_INLINE_FUNCTION static void store(View out,std::size_t t,const double* sum)
        { int j=0; for(int a=0;a<3;++a) {
      for(int b=a;b<3;++b) { out(t,3*a+b)=out(t,3*b+a)=sum[j++]; }
    }
  }
    };
    template<class ExecutionSpace=Kokkos::DefaultExecutionSpace>
    class KokkosStressFMM : public KokkosParticleFMM<9,6,ExecutionSpace,FMMStressMatrixOutput>
    {
        using Base=KokkosParticleFMM<9,6,ExecutionSpace,FMMStressMatrixOutput>;
    public:
        using Base::Base;
        template<class Kernel,class PositionsView,class TensorView,class OutputView>
        void evaluate(const ExecutionSpace& ex,const Kernel& kernel,PositionsView src,TensorView tensor,
                      PositionsView tgt,OutputView output,std::size_t ns,std::size_t nt)
        { Base::evaluate(ex,FMMStressKernelAdapter<Kernel>(kernel),src,tensor,tgt,output,ns,nt); }
    };
}
#endif
