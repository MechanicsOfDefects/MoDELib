/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosShortRangeParticleNeighbors_H_
#define model_KokkosShortRangeParticleNeighbors_H_
#include <cfloat>
#include <Kokkos_Core.hpp>
#include <FMM/ShortRangeParticleNeighbors.h>
namespace model
{
    template<class ExecutionSpace=Kokkos::DefaultExecutionSpace>
    class KokkosShortRangeParticleNeighbors
    {
        using Indices=Kokkos::View<std::size_t*,ExecutionSpace>;
        using Shifts=Kokkos::View<double**,Kokkos::LayoutLeft,ExecutionSpace>;
        Indices offsets_,sources_,images_;
        Shifts shifts_;
        std::size_t ns_,nt_;
        double cutoff2_;
    public:
        // Own the uploaded list. Particle buffers must still refer to the same snapshot.
        KokkosShortRangeParticleNeighbors(const ExecutionSpace& ex,const ShortRangeParticleNeighbors& list)
            :offsets_("neighbor offsets",list.offsets.size()),sources_("neighbor sources",list.pairs.size()),
             images_("neighbor images",list.pairs.size()),shifts_("neighbor shifts",list.shifts.size(),3),
             ns_(list.sourceCount),nt_(list.targetCount),cutoff2_(list.cutoff*list.cutoff)
        {
            auto o=Kokkos::create_mirror(offsets_),s=Kokkos::create_mirror(sources_),im=Kokkos::create_mirror(images_);
            auto sh=Kokkos::create_mirror(shifts_);
            for(std::size_t i=0;i<list.offsets.size();++i) { o(i)=list.offsets[i]; }
            for(std::size_t i=0;i<list.pairs.size();++i) { s(i)=list.pairs[i].source; im(i)=list.pairs[i].image; }
            for(std::size_t i=0;i<list.shifts.size();++i) { for(int j=0;j<3;++j) { sh(i,j)=list.shifts[i][j]; } }
            Kokkos::deep_copy(ex,offsets_,o); Kokkos::deep_copy(ex,sources_,s);
            Kokkos::deep_copy(ex,images_,im); Kokkos::deep_copy(ex,shifts_,sh); ex.fence();
        }
        template<class Kernel,class Positions,class Tensors,class Output>
        void evaluate(const ExecutionSpace& ex,const Kernel& kernel,Positions src,Tensors tensor,Positions tgt,Output stress,
                      std::size_t sourceCount,std::size_t targetCount) const
        {
            if(sourceCount!=ns_||targetCount!=nt_||src.extent(0)<ns_||tgt.extent(0)<nt_||tensor.extent(0)<ns_||stress.extent(0)<nt_||
               src.extent(1)!=3||tgt.extent(1)!=3||tensor.extent(1)!=9||stress.extent(1)!=9)
            { throw std::invalid_argument("Neighbor list / particle snapshot shape mismatch"); }
            const auto offsets=offsets_,sources=sources_,images=images_; const auto shifts=shifts_;
            const double cutoff2=cutoff2_;
            Kokkos::parallel_for("Short range stress neighbors",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex,0,nt_),
                KOKKOS_LAMBDA(const std::size_t t)
                {
                    double sum[6]={};
                    for(std::size_t p=offsets(t);p<offsets(t+1);++p)
                    {
                        const auto s=sources(p),im=images(p); double r[3],q[9],d2=0;
                        for(int j=0;j<3;++j) { r[j]=(tgt(t,j)-src(s,j))-shifts(im,j); d2+=r[j]*r[j]; }
                        if(d2<=cutoff2)
                        {
                            for(int j=0;j<9;++j) { q[j]=tensor(s,j); }
                            kernel.addSymmetric(r,q,sum);
                        }
                    }
                    int p=0;
                    for(int i=0;i<3;++i) { for(int j=i;j<3;++j) { stress(t,3*i+j)=stress(t,3*j+i)=sum[p++]; } }
                });
        }
    };
}
#endif
