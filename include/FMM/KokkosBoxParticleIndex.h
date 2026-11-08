/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosBoxParticleIndex_H_
#define model_KokkosBoxParticleIndex_H_
#include <cfloat>
#include <Kokkos_Core.hpp>
#include <FMM/BoxParticleIndex.h>
namespace model
{
    // Upload after rebuilding the host index. Views contain original point IDs.
    template<int dim, class ExecutionSpace=Kokkos::DefaultExecutionSpace>
    class KokkosBoxParticleIndex
    {
        using View=Kokkos::View<std::size_t*,ExecutionSpace>;
        View sourceOffsets_,sourcePermutation_,targetOffsets_,targetPermutation_;
        static View copy(const std::vector<std::size_t>& input,const std::string& label)
        {
            View result(label,input.size());
            auto host=Kokkos::create_mirror(result);
            for(std::size_t i=0;i<input.size();++i)
            {
                host(i)=input[i];
            }
            Kokkos::deep_copy(result,host);
            return result;
        }
    public:
        using ConstView=Kokkos::View<const std::size_t*,ExecutionSpace>;
        ConstView sourceOffsets() const { return sourceOffsets_; }
        ConstView sourcePermutation() const { return sourcePermutation_; }
        ConstView targetOffsets() const { return targetOffsets_; }
        ConstView targetPermutation() const { return targetPermutation_; }
        void upload(const BoxParticleIndex<dim>& index,const BoxTree<dim>& tree)
        {
            if(!Kokkos::is_initialized() || Kokkos::is_finalized())
            {
                throw std::logic_error("Initialize Kokkos before uploading leaf indices");
            }
            if(!index.matches(tree))
            {
                throw std::logic_error("Stale BoxParticleIndex");
            }
            auto so=copy(index.sources().offsets,"Source leaf offsets");
            auto sp=copy(index.sources().permutation,"Source leaf point IDs");
            auto to=copy(index.targets().offsets,"Target leaf offsets");
            auto tp=copy(index.targets().permutation,"Target leaf point IDs");
            sourceOffsets_=so; sourcePermutation_=sp;
            targetOffsets_=to; targetPermutation_=tp;
        }
    };
}
#endif
