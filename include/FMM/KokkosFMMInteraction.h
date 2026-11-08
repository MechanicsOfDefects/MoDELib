/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosFMMInteraction_H_
#define model_KokkosFMMInteraction_H_
#include <FMM/KokkosParticleFMM.h>
#include <FMM/FMMManager.h>
#include <FMM/ShortRangeParticleNeighbors.h>
#include <FMM/LeafBoxNeighbors.h>
#include <memory>
namespace model
{
    enum class FMMBoundaryMode { freeSpace, periodic };
    struct FMMInteractionOptions
    {
        int interpolationOrder=4;
        bool boxShortRange=true; // Shared leaf stencil; false only for legacy reference tests.
        bool screened=true;
        bool interpolateScreenedNear=true; // False is a direct-reference/debug path.
        bool externalNear=false; // Caller supplies full segment near interactions in unscreened mode.
        bool buildShortNeighbors=true; // False for callers assembling their own sparse local operator.
        double shortCutoff=1; // Legacy radial reference path only; ignored with boxShortRange.
        std::size_t translationCacheBytes=256u*1024u*1024u;
        FMMBoundaryMode boundary=FMMBoundaryMode::freeSpace;
    };
    // Geometry snapshots are owned, never borrowed from mutable physics arrays.
    // Source strengths and outputs stay caller-owned. No concurrent evaluations on one plan.
    template<class Family,class ExecutionSpace=Kokkos::DefaultExecutionSpace>
    class KokkosFMMInteraction
    {
    public:
        using Positions=std::array<std::vector<double>,3>;
        using View=Kokkos::View<double**,Kokkos::LayoutLeft,ExecutionSpace>;
        static constexpr int NS=Family::sourceComponents,NF=Family::fieldComponents;
    private:
        using Plan=KokkosParticleFMM<NS,NF,ExecutionSpace>;
        using Indices=Kokkos::View<std::size_t*,ExecutionSpace>;
        struct State
        {
            std::shared_ptr<const BoxTree<3>> tree;
            std::uint64_t revision;
            std::size_t ns,nt,shortPairs=0;
            View sources,targets,scratch;
            Indices offsets,ids,sourceOffsets,sourceIDs,targetLeaves;
            Plan plan;
            static const BoxTree<3>& checkedTree(const std::shared_ptr<const BoxTree<3>>& t)
            { if(!t) { throw std::invalid_argument("Null FMM tree"); } return *t; }
            State(const ExecutionSpace& ex,std::shared_ptr<const BoxTree<3>> t,const Positions& src,const Positions& tgt,
                  const Family& family,const FMMInteractionOptions& opt)
                :tree(std::move(t)),revision(checkedTree(tree).revision()),ns(src[0].size()),nt(tgt[0].size()),
                 sources("FMM source snapshot",ns,3),targets("FMM target snapshot",nt,3),scratch("FMM short result",nt,NF),
                 plan(ex,*tree,src,tgt,opt.interpolationOrder)
            {
                auto sh=Kokkos::create_mirror(sources),th=Kokkos::create_mirror(targets);
                for(int c=0;c<3;++c)
                {
                    for(std::size_t i=0;i<ns;++i) { sh(i,c)=src[c][i]; }
                    for(std::size_t i=0;i<nt;++i) { th(i,c)=tgt[c][i]; }
                }
                Kokkos::deep_copy(ex,sources,sh); Kokkos::deep_copy(ex,targets,th);
                if(opt.screened && opt.buildShortNeighbors)
                {
                if(opt.boxShortRange)
                {
                    LeafBoxNeighbors list(*tree,src,tgt);
                    const auto upload=[&](const char* name,const std::vector<std::size_t>& values){Indices v(std::string(name),values.size());auto h=Kokkos::create_mirror(v);for(std::size_t i=0;i<values.size();++i) {
              h(i)=values[i];
            }
            Kokkos::deep_copy(ex,v,h);return v;};
                    for(std::size_t b=0;b<list.offsets.size()-1;++b)
                    { const auto targetsInBox=list.index.targets().offsets[b+1]-list.index.targets().offsets[b];
                      for(auto k=list.offsets[b];k<list.offsets[b+1];++k) {
              shortPairs+=targetsInBox*(list.index.sources().offsets[list.boxes[k]+1]-list.index.sources().offsets[list.boxes[k]]);
            }
          }
                    offsets=upload("short box offsets",list.offsets);ids=upload("short boxes",list.boxes);
                    sourceOffsets=upload("short source offsets",list.index.sources().offsets);
                    sourceIDs=upload("short source IDs",list.index.sources().permutation);
                    targetLeaves=upload("short target leaves",list.index.targets().leafOfPoint);
                }
                else {
                ShortRangeParticleNeighbors list(src,tgt,family.searchRadius(opt.shortCutoff));
                shortPairs=list.pairs.size();
                offsets=Indices("FMM short offsets",list.offsets.size()); ids=Indices("FMM short sources",list.pairs.size());
                auto oh=Kokkos::create_mirror(offsets),ih=Kokkos::create_mirror(ids);
                for(std::size_t i=0;i<list.offsets.size();++i) { oh(i)=list.offsets[i]; }
                for(std::size_t i=0;i<list.pairs.size();++i) { ih(i)=list.pairs[i].source; }
                Kokkos::deep_copy(ex,offsets,oh); Kokkos::deep_copy(ex,ids,ih);
                }
                }
                plan.setExternalNear(opt.externalNear);
                plan.setInterpolateNearBoxes(opt.screened && opt.interpolateScreenedNear);
                plan.setTranslationCacheBudget(opt.translationCacheBytes); ex.fence();
            }
        };
        std::shared_ptr<const Family> family_;
        ExecutionSpace ex_;
        FMMInteractionOptions options_;
        std::unique_ptr<State> state_;
        template<class Strengths,class Output> void check(Strengths q,Output out) const
        {
            if(state_->tree->revision()!=state_->revision) { throw std::logic_error("FMM geometry changed; rebuildGeometry required"); }
            if(q.extent(0)!=state_->ns || q.extent(1)!=NS || out.extent(0)!=state_->nt || out.extent(1)!=NF)
            { throw std::invalid_argument("FMM interaction source/output shape mismatch"); }
            if(static_cast<const void*>(q.data())==static_cast<const void*>(out.data()) && q.size())
            { throw std::invalid_argument("FMM source and output must not alias"); }
        }
    public:
        KokkosFMMInteraction(std::shared_ptr<const Family> family,const ExecutionSpace& ex,std::shared_ptr<const BoxTree<3>> tree,
                            const Positions& src,const Positions& tgt,FMMInteractionOptions options={})
            :family_(std::move(family)),ex_(ex),options_(options)
        {
            if(!family_) { throw std::invalid_argument("Null FMM kernel family"); }
            if(options_.boundary!=FMMBoundaryMode::freeSpace)
            { throw std::invalid_argument("Periodic FMM correction is not implemented; refusing free-space substitution"); }
            if(options_.screened && !options_.boxShortRange && (!std::isfinite(options_.shortCutoff)||options_.shortCutoff<=0))
            { throw std::invalid_argument("FMM short cutoff must be finite and positive"); }
            rebuildGeometry(std::move(tree),src,tgt);
        }
        void rebuildGeometry(std::shared_ptr<const BoxTree<3>> tree,const Positions& src,const Positions& tgt)
        {
            ex_.fence();
            auto replacement=std::make_unique<State>(ex_,std::move(tree),src,tgt,*family_,options_);
            state_=std::move(replacement);
        }
        template<class Strengths,class Output> void evaluateLong(Strengths strengths,Output out)
        {
            check(strengths,out);
            if(!options_.screened)
            { state_->plan.evaluate(ex_,family_->fullKernel(),state_->sources,strengths,state_->targets,out,state_->ns,state_->nt); return; }
            state_->plan.evaluate(ex_,family_->longKernel(),state_->sources,strengths,state_->targets,out,state_->ns,state_->nt);
        }
        template<class Strengths,class Output> void evaluateShort(Strengths strengths,Output out)
        {
            if(!options_.screened) {
      throw std::logic_error("Unscreened FMM has no additive short-range correction");
    }
    if(!options_.buildShortNeighbors) {
      throw std::logic_error("This interaction was prepared for long-range evaluation only");
    }
    check(strengths,out);
            const auto src=state_->sources,tgt=state_->targets;
            const auto offsets=state_->offsets,ids=state_->ids;
            const auto kernel=family_->shortKernel(); const double cutoff=options_.shortCutoff;
            const bool boxShort=options_.boxShortRange;
            const auto so=state_->sourceOffsets,si=state_->sourceIDs,tl=state_->targetLeaves;
            Kokkos::parallel_for("FMM registered short kernel",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex_,0,state_->nt),KOKKOS_LAMBDA(const std::size_t t)
            {
                double sum[NF]={};
                const auto b=boxShort?tl(t):t;
                for(std::size_t j=offsets(b);j<offsets(b+1);++j)
                {
                  const auto begin=boxShort?so(ids(j)):j,end=boxShort?so(ids(j)+1):j+1;
                  for(auto a=begin;a<end;++a){
                    const auto s=boxShort?si(a):ids(a); double r[3],q[NS];
                    for(int c=0;c<3;++c) { r[c]=tgt(t,c)-src(s,c); }
                    if(!boxShort && kernel.distance(r)>cutoff) { continue; }
                    for(int c=0;c<NS;++c) { q[c]=strengths(s,c); }
                    kernel.add(r,q,sum);
                  }
                }
                for(int c=0;c<NF;++c) { out(t,c)=sum[c]; }
            });
            ex_.fence();
        }
        template<class Strengths,class Output> void evaluate(Strengths strengths,Output out)
        {
            evaluateLong(strengths,out);
            if(!options_.screened) {
      return;
    }
    evaluateShort(strengths,state_->scratch);
            const auto local=state_->scratch;
            Kokkos::parallel_for("FMM combine split",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex_,0,state_->nt),KOKKOS_LAMBDA(const std::size_t t)
            { for(int c=0;c<NF;++c) { out(t,c)+=local(t,c); } });
            ex_.fence();
        }
        void setIntegratedLeafMoments(View values)
        {state_->plan.setIntegratedLeafMoments(ex_,values);}
        std::vector<double> hostLeafLocals() const { return state_->plan.hostLeafLocals(); }
        std::size_t translationCacheBytes() const { return state_->plan.translationCacheBytes(); }
        std::size_t shortPairCount() const { return state_->shortPairs; }
    };
}
#endif
