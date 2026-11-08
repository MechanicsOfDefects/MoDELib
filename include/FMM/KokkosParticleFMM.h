/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosParticleFMM_H_
#define model_KokkosParticleFMM_H_
#include <cfloat>
#include <Kokkos_Core.hpp>
#include <FMM/BoxParticleIndex.h>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <chrono>
#include <type_traits>
#include <utility>
#include <typeindex>
#include <limits>
namespace model
{
    template<int N> struct FMMCompactOutput
    {
        static constexpr int components=N;
        template<class View> KOKKOS_INLINE_FUNCTION static void store(View output,std::size_t t,const double* sum)
        { for(int c=0;c<N;++c) { output(t,c)=sum[c]; } }
    };
    template<class Kernel,class=void> struct HasFMMTranslationCache:std::false_type {};
    template<class Kernel> struct HasFMMTranslationCache<Kernel,std::void_t<decltype(std::declval<Kernel>().translationKey())>>:std::true_type {};
    // Free-space interpolation FMM for arbitrary fixed-size linear particle kernels.
    // Geometry/weights are a snapshot: rebuild after point motion/reordering.
    // P2M, M2M, M2L, L2L, L2P and exact long-range near interactions.
    template<int SourceComponents,int FieldComponents,class ExecutionSpace=Kokkos::DefaultExecutionSpace,
             class OutputPolicy=FMMCompactOutput<FieldComponents>>
    class KokkosParticleFMM
    {
        static_assert(SourceComponents>0 && FieldComponents>0,"Positive FMM component counts required");
        using View=Kokkos::View<double**,Kokkos::LayoutLeft,ExecutionSpace>;
        using Indices=Kokkos::View<std::size_t*,ExecutionSpace>;
        using Positions=std::array<std::vector<double>,3>;
        const BoxTree<3>* tree_;
        std::uint64_t revision_;
        unsigned depth_;
        int order_,nodes_;
        std::size_t ns_,nt_,boxCount_,leafStart_;
        std::vector<std::size_t> levelStart_;
        View centers_,nodeCoordinates_,sourceWeights_,targetWeights_,transfer_,moments_,locals_;
        Indices sourceOffsets_,sourceIDs_,targetLeaves_,m2lOffsets_,m2lSources_,nearOffsets_,nearSources_,sourceCounts_,targetCounts_;
        std::size_t farPairs_=0,nearPairs_=0;
        bool externalNear_=false,interpolateNear_=false,interpolateNearBoxes_=false,profile_=false;
        std::array<double,4> lastTimings_{}; // cache preparation, upward, downward, target/near (ms)
        std::size_t cacheBudget_=256u*1024u*1024u;
        View translations_;
        View integratedLeafMoments_;
        Indices translationIDs_,translationTargets_,translationSources_;
        std::vector<std::size_t> hostTranslationIDs_,hostTranslationTargets_,hostTranslationSources_;
        std::vector<double> translationKey_{};
        std::type_index translationType_{typeid(void)};
        bool cacheReady_=false;
        template<class Kernel> void prepareTranslations(const ExecutionSpace& ex,const Kernel& kernel)
        {
            if constexpr(HasFMMTranslationCache<Kernel>::value)
            {
                const auto rawKey=kernel.translationKey();
                const std::vector<double> key(rawKey.begin(),rawKey.end());
                if(cacheReady_ && translationType_==std::type_index(typeid(Kernel)) && translationKey_==key) { return; }
                ex.fence();
                const std::size_t perClass=std::size_t(nodes_)*nodes_;
                const auto count=std::min(hostTranslationTargets_.size(),cacheBudget_/(perClass*Kernel::translationComponents*sizeof(double)));
                auto ids=hostTranslationIDs_;
                for(auto& id:ids) { if(id>=count) { id=std::numeric_limits<std::size_t>::max(); } }
                translationIDs_=upload(ex,"FMM translation IDs",ids);
                translations_=View("FMM cached translation coefficients",count*perClass,Kernel::translationComponents);
                const auto xyz=nodeCoordinates_,cache=translations_;
                const auto targets=translationTargets_,sources=translationSources_; const int P=nodes_;
                Kokkos::parallel_for("FMM prepare translations",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex,0,count*perClass),KOKKOS_LAMBDA(const std::size_t k)
                {
                    const auto id=k/(P*P),j=(k/P)%P,a=k%P;
                    double r[3],data[Kernel::translationComponents];
                    for(int c=0;c<3;++c) { r[c]=xyz(targets(id)*P+j,c)-xyz(sources(id)*P+a,c); }
                    kernel.prepareTranslation(r,data);
                    for(int c=0;c<Kernel::translationComponents;++c) { cache(k,c)=data[c]; }
                });
                ex.fence(); translationKey_=key; translationType_=std::type_index(typeid(Kernel)); cacheReady_=true;
            }
        }
        static int checkedOrder(int order)
        {
            if(order<2||order>6) { throw std::invalid_argument("FMM interpolation order must be 2..6"); }
            return order;
        }
        static std::vector<double> weights(double x,const std::vector<double>& nodes)
        {
            std::vector<double> w(nodes.size());
            for(std::size_t i=0;i<nodes.size();++i)
            {
                if(x==nodes[i]) { w[i]=1; return w; }
            }
            double sum=0;
            for(std::size_t i=0;i<nodes.size();++i)
            {
                w[i]=(i%2?-1.0:1.0)*((i==0||i+1==nodes.size())?.5:1)/(x-nodes[i]); sum+=w[i];
            }
            for(auto& v:w) { v/=sum; }
            return w;
        }
        static Indices upload(const ExecutionSpace& ex,const char* label,const std::vector<std::size_t>& a)
        {
            Indices v(std::string(label),a.size()); auto h=Kokkos::create_mirror(v);
            for(std::size_t i=0;i<a.size();++i) { h(i)=a[i]; }
            Kokkos::deep_copy(ex,v,h); ex.fence(); return v;
        }
    public:
        KokkosParticleFMM(const ExecutionSpace& ex,const BoxTree<3>& tree,const Positions& sources,const Positions& targets,int order=4)
            :tree_(&tree),revision_(tree.revision()),depth_(tree.depth()),order_(checkedOrder(order)),nodes_(order_*order_*order_),ns_(sources[0].size()),nt_(targets[0].size())
        {
            BoxParticleIndex<3> index; index.rebuild(tree,sources,targets);
            boxCount_=0;
            for(unsigned l=0;l<=depth_;++l) { levelStart_.push_back(boxCount_); boxCount_+=std::size_t(1)<<(3*l); }
            leafStart_=levelStart_.back(); const std::size_t leaves=std::size_t(1)<<(3*depth_);
            std::vector<double> nodes(order);
            for(int i=0;i<order;++i) { nodes[i]=std::cos(3.14159265358979323846*i/(order-1)); }
            centers_=View("FMM centers",boxCount_,3); nodeCoordinates_=View("FMM interpolation positions",boxCount_*nodes_,3);
            auto ch=Kokkos::create_mirror(centers_),nh=Kokkos::create_mirror(nodeCoordinates_);
            std::vector<Eigen::Vector3d> center(boxCount_); std::vector<double> diameter(depth_+1);
            for(unsigned l=0;l<=depth_;++l)
            {
                const int n=1<<l; const Eigen::Vector3d width=(tree.upper()-tree.lower())/n; diameter[l]=width.norm();
                for(int z=0;z<n;++z) { for(int y=0;y<n;++y) { for(int x=0;x<n;++x)
                {
                    const auto b=levelStart_[l]+x+n*(y+n*z); center[b]=tree.lower()+width.cwiseProduct(Eigen::Vector3d(x+.5,y+.5,z+.5));
                    for(int a=0;a<3;++a) { ch(b,a)=center[b][a]; }
                    for(int j=0;j<nodes_;++j)
                    {
                        const int ij[3]={j%order,(j/order)%order,j/(order*order)};
                        for(int a=0;a<3;++a) { nh(b*nodes_+j,a)=center[b][a]+.5*width[a]*nodes[ij[a]]; }
                    }
                } } }
            }
            Kokkos::deep_copy(ex,centers_,ch); Kokkos::deep_copy(ex,nodeCoordinates_,nh);
            std::vector<std::size_t> boxToLex(leaves);
            for(std::size_t i=0;i<leaves;++i)
            {
                const auto key=tree.boxes(depth_)[i]->key(); const auto cellWidth=BoxTree<3>::latticeSize/(1u<<depth_); const auto n=1u<<depth_;
                boxToLex[i]=key[0]/cellWidth+n*(key[1]/cellWidth+n*(key[2]/cellWidth));
            }
            std::vector<std::vector<std::size_t>> grouped(leaves);
            std::vector<std::size_t> sc(boxCount_,0),tc(boxCount_,0),tl(nt_);
            for(std::size_t s=0;s<ns_;++s) { const auto b=boxToLex[index.sources().leafOfPoint[s]]; grouped[b].push_back(s); ++sc[leafStart_+b]; }
            for(std::size_t t=0;t<nt_;++t) { tl[t]=boxToLex[index.targets().leafOfPoint[t]]; ++tc[leafStart_+tl[t]]; }
            for(unsigned l=depth_;l>0;--l)
            {
                const int n=1<<l,p=n/2;
                for(int z=0;z<n;++z) { for(int y=0;y<n;++y) { for(int x=0;x<n;++x)
                {
                    const auto b=levelStart_[l]+x+n*(y+n*z),parent=levelStart_[l-1]+x/2+p*(y/2+p*(z/2));
                    sc[parent]+=sc[b]; tc[parent]+=tc[b];
                } } }
            }
            sourceCounts_=upload(ex,"FMM source counts",sc); targetCounts_=upload(ex,"FMM target counts",tc);
            targetLeaves_=upload(ex,"FMM target leaves",tl);
            std::vector<std::size_t> so(leaves+1),si;
            for(std::size_t b=0;b<leaves;++b) { si.insert(si.end(),grouped[b].begin(),grouped[b].end()); so[b+1]=si.size(); }
            sourceOffsets_=upload(ex,"FMM source offsets",so); sourceIDs_=upload(ex,"FMM source IDs",si);
            sourceWeights_=View("FMM source weights",ns_,3*order); targetWeights_=View("FMM target weights",nt_,3*order);
            const Eigen::Vector3d leafWidth=(tree.upper()-tree.lower())/(1u<<depth_);
            const auto setWeights=[&](const Positions& pts,View view,const std::vector<std::size_t>& leafOfPoint)
            {
                auto host=Kokkos::create_mirror(view);
                for(std::size_t p=0;p<pts[0].size();++p) { for(int a=0;a<3;++a)
                {
                    const auto w=weights(2*(pts[a][p]-center[leafStart_+leafOfPoint[p]][a])/leafWidth[a],nodes);
                    for(int j=0;j<order;++j) { host(p,a*order+j)=w[j]; }
                } }
                Kokkos::deep_copy(ex,view,host); ex.fence();
            };
            std::vector<std::size_t> sl(ns_); for(std::size_t i=0;i<ns_;++i) { sl[i]=boxToLex[index.sources().leafOfPoint[i]]; }
            setWeights(sources,sourceWeights_,sl); setWeights(targets,targetWeights_,tl);
            // Parent cardinal basis evaluated at child nodes, for each octant.
            transfer_=View("FMM child transfer",8*nodes_,nodes_); auto th=Kokkos::create_mirror(transfer_);
            for(int oct=0;oct<8;++oct) { for(int child=0;child<nodes_;++child)
            {
                const int ci[3]={child%order,(child/order)%order,child/(order*order)};
                std::vector<double> w[3];
                for(int a=0;a<3;++a) { w[a]=weights(.5*nodes[ci[a]]+((oct>>a)&1?.5:-.5),nodes); }
                for(int parent=0;parent<nodes_;++parent) { th(oct*nodes_+parent,child)=w[0][parent%order]*w[1][(parent/order)%order]*w[2][parent/(order*order)]; }
            } }
            Kokkos::deep_copy(ex,transfer_,th);
            // Dual traversal. Each source-target box pair is accepted once or refined.
            // Sphere separation also handles highly rectangular root boxes safely.
            std::vector<std::vector<std::size_t>> far(boxCount_),near(leaves);
            std::function<void(unsigned,std::size_t,std::size_t)> walk;
            const auto childID=[&](unsigned l,std::size_t b,int oct)
            {
                const std::size_t n=1u<<l,m=2*n,local=b-levelStart_[l];
                const auto x=local%n,y=(local/n)%n,z=local/(n*n);
                return levelStart_[l+1]+2*x+(oct&1)+m*(2*y+((oct>>1)&1)+m*(2*z+((oct>>2)&1)));
            };
            walk=[&](unsigned l,std::size_t t,std::size_t s)
            {
                if(!tc[t]||!sc[s]) { return; }
                if((center[t]-center[s]).norm()>1.01*diameter[l]) { far[t].push_back(s); ++farPairs_; return; }
                if(l==depth_) { near[t-leafStart_].push_back(s-leafStart_); ++nearPairs_; return; }
                for(int a=0;a<8;++a) { for(int b=0;b<8;++b) { walk(l+1,childID(l,t,a),childID(l,s,b)); } }
            };
            walk(0,0,0);
            const auto pack=[&](const auto& rows,const char* ol,const char* il,Indices& offsets,Indices& ids)
            {
                std::vector<std::size_t> o(rows.size()+1),v;
                for(std::size_t b=0;b<rows.size();++b) { v.insert(v.end(),rows[b].begin(),rows[b].end()); o[b+1]=v.size(); }
                offsets=upload(ex,ol,o); ids=upload(ex,il,v);
            };
            pack(far,"FMM M2L offsets","FMM M2L sources",m2lOffsets_,m2lSources_);
            pack(near,"FMM near offsets","FMM near sources",nearOffsets_,nearSources_);
            // Translation invariance: share one operator for each level and integer displacement.
            std::map<std::array<int,4>,std::size_t> classes;
            for(unsigned l=0;l<=depth_;++l)
            {
                const int n=1<<l;
                for(std::size_t b=levelStart_[l];b<levelStart_[l]+std::size_t(n)*n*n;++b)
                {
                    const int t=static_cast<int>(b-levelStart_[l]);
                    for(const auto source:far[b])
                    {
                        const int s=static_cast<int>(source-levelStart_[l]);
                        const std::array<int,4> key{static_cast<int>(l),t%n-s%n,(t/n)%n-(s/n)%n,t/(n*n)-s/(n*n)};
                        const auto result=classes.emplace(key,classes.size());
                        if(result.second) { hostTranslationTargets_.push_back(b); hostTranslationSources_.push_back(source); }
                        hostTranslationIDs_.push_back(result.first->second);
                    }
                }
            }
            translationTargets_=upload(ex,"FMM translation targets",hostTranslationTargets_);
            translationSources_=upload(ex,"FMM translation sources",hostTranslationSources_);
            moments_=View("FMM multipoles",boxCount_*nodes_,SourceComponents); locals_=View("FMM locals",boxCount_*nodes_,FieldComponents);
            ex.fence();
        }
        KokkosParticleFMM(const KokkosParticleFMM&)=delete;
        KokkosParticleFMM& operator=(const KokkosParticleFMM&)=delete;
        // Snapshot of leaf local coefficients after evaluate(); safe for read-only host queries.
        std::vector<double> hostLeafLocals() const
        {
            auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),locals_);
            std::vector<double> values((boxCount_-leafStart_)*nodes_*FieldComponents);
            for(std::size_t b=leafStart_;b<boxCount_;++b) {
      for(int j=0;j<nodes_;++j) {
        for(int c=0;c<FieldComponents;++c) {
          values[((b-leafStart_)*nodes_+j)*FieldComponents+c]=host(b*nodes_+j,c);
        }
      }
    }
    return values;
        }
        // Coefficients use lexical leaf order, then tensor-product node order.
        // The geometry snapshot must contain a marker in every occupied source leaf.
        void setIntegratedLeafMoments(const ExecutionSpace& ex,View values)
        {
            if(!interpolateNearBoxes_&&!externalNear_) {
      throw std::logic_error("Integrated sources require interpolated near boxes");
    }
    if(values.extent(0)!=(boxCount_-leafStart_)*nodes_||values.extent(1)!=SourceComponents) {
      throw std::invalid_argument("Integrated leaf coefficient shape mismatch");
    }
    if(integratedLeafMoments_.extent(0)!=values.extent(0)||integratedLeafMoments_.extent(1)!=SourceComponents)
    {
        integratedLeafMoments_=View("integrated source coefficients",values.extent(0),SourceComponents);
    }
            Kokkos::deep_copy(ex,integratedLeafMoments_,values);ex.fence();
        }
        std::size_t farBoxPairs() const { return farPairs_; }
        std::size_t nearBoxPairs() const { return nearPairs_; }
        // Near interpolation is experimental and has no automatic error bound.
        // Validate against direct evaluation for each geometry and kernel scale.
        void setProfiling(bool enabled) { profile_=enabled; }
        const std::array<double,4>& lastTimings() const { return lastTimings_; }
        void setExternalNear(bool enabled){externalNear_=enabled;}
        void setInterpolateNearBoxes(bool enabled) { interpolateNearBoxes_=enabled; }
        void setInterpolateNear(bool enabled) { interpolateNear_=enabled; }
        void setTranslationCacheBudget(std::size_t bytes) { cacheBudget_=bytes; cacheReady_=false; }
        std::size_t translationCacheBytes() const { return translations_.size()*sizeof(double); }
        // Geometry is fixed, source tensor strengths may change between evaluations.
        template<class Kernel,class PositionsView,class TensorView,class OutputView>
        void evaluate(const ExecutionSpace& ex,const Kernel& kernel,PositionsView src,TensorView tensor,PositionsView tgt,OutputView output,
                      std::size_t ns,std::size_t nt)
        {
            if(tree_->revision()!=revision_) { throw std::logic_error("FMM tree changed; rebuild the plan"); }
            if(ns!=ns_||nt!=nt_||src.extent(0)<ns||tensor.extent(0)<ns||tgt.extent(0)<nt||output.extent(0)<nt||
               (ns && (src.extent(1)!=3||tensor.extent(1)!=SourceComponents))||(nt && (tgt.extent(1)!=3||output.extent(1)!=OutputPolicy::components)))
            { throw std::invalid_argument("FMM particle snapshot shape mismatch"); }
            using Clock=std::chrono::steady_clock;
            auto stamp=Clock::now();
            const auto checkpoint=[&](int stage)
            {
                if(profile_) { ex.fence(); const auto now=Clock::now(); lastTimings_[stage]=std::chrono::duration<double,std::milli>(now-stamp).count(); stamp=now; }
            };
            prepareTranslations(ex,kernel); checkpoint(0);
            const bool externalNear=externalNear_,interpolateNear=interpolateNear_,interpolateNearBoxes=interpolateNearBoxes_;
            const auto cache=translations_; const auto cacheIDs=translationIDs_;
            const int p=order_,P=nodes_; const auto leaf=leafStart_;
            const auto M=moments_,L=locals_,sw=sourceWeights_,tw=targetWeights_,tr=transfer_,xyz=nodeCoordinates_;
            const auto so=sourceOffsets_,si=sourceIDs_,tl=targetLeaves_,sc=sourceCounts_,tc=targetCounts_;
            const auto fo=m2lOffsets_,fi=m2lSources_,no=nearOffsets_,ni=nearSources_;
            Kokkos::deep_copy(ex,M,0.0); Kokkos::deep_copy(ex,L,0.0);
            const auto leaves=std::size_t(1)<<(3*depth_);
            if(integratedLeafMoments_.size())
            {
                const auto integrated=integratedLeafMoments_;
                Kokkos::parallel_for("FMM integrated segment moments",Kokkos::RangePolicy<ExecutionSpace>(ex,0,leaves*P),KOKKOS_LAMBDA(const int k)
                {for(int c=0;c<SourceComponents;++c) {
              M(leaf*P+k,c)=integrated(k,c);
            }
          });
            }
            else {
      Kokkos::parallel_for("FMM P2M",Kokkos::RangePolicy<ExecutionSpace>(ex,0,leaves*P),KOKKOS_LAMBDA(const int k)
            {
                const auto b=k/P; const int j=k%P; double sum[SourceComponents]={};
                for(std::size_t a=so(b);a<so(b+1);++a)
                {
                    const auto s=si(a); const double w=sw(s,j%p)*sw(s,p+(j/p)%p)*sw(s,2*p+j/(p*p));
                    for(int c=0;c<SourceComponents;++c) { sum[c]+=w*tensor(s,c); }
                }
                for(int c=0;c<SourceComponents;++c) { M((leaf+b)*P+j,c)=sum[c]; }
            });
    }
    for(unsigned l=depth_;l>0;--l)
            {
                const auto parentStart=levelStart_[l-1],childStart=levelStart_[l]; const int n=1<<(l-1),m=2*n;
                Kokkos::parallel_for("FMM M2M",Kokkos::RangePolicy<ExecutionSpace>(ex,0,n*n*n*P),KOKKOS_LAMBDA(const int k)
                {
                    const int b=k/P,j=k%P,x=b%n,y=(b/n)%n,z=b/(n*n); if(!sc(parentStart+b)) { return; }
                    double sum[SourceComponents]={};
                    for(int oct=0;oct<8;++oct)
                    {
                        const auto child=childStart+2*x+(oct&1)+m*(2*y+((oct>>1)&1)+m*(2*z+((oct>>2)&1)));
                        if(!sc(child)) { continue; }
                        for(int a=0;a<P;++a) { const double w=tr(oct*P+j,a); for(int c=0;c<SourceComponents;++c) { sum[c]+=w*M(child*P+a,c); } }
                    }
                    for(int c=0;c<SourceComponents;++c) { M((parentStart+b)*P+j,c)=sum[c]; }
                });
            }
            checkpoint(1);
            for(unsigned l=1;l<=depth_;++l)
            {
                const auto start=levelStart_[l]; const int n=1<<l;
                Kokkos::parallel_for("FMM M2L",Kokkos::RangePolicy<ExecutionSpace>(ex,0,n*n*n*P),KOKKOS_LAMBDA(const int k)
                {
                    const auto b=start+k/P; const int j=k%P; if(!tc(b)) { return; }
                    double sum[FieldComponents]={};
                    for(std::size_t pair=fo(b);pair<fo(b+1);++pair)
                    {
                        const auto s=fi(pair);
                        for(int a=0;a<P;++a)
                        {
                            double r[3],q[SourceComponents]; for(int c=0;c<3;++c) { r[c]=xyz(b*P+j,c)-xyz(s*P+a,c); }
                            for(int c=0;c<SourceComponents;++c) { q[c]=M(s*P+a,c); }
                            if constexpr(HasFMMTranslationCache<Kernel>::value)
                            {
                                const auto id=cacheIDs(pair);
                                if(id!=std::numeric_limits<std::size_t>::max())
                                {
                                    double data[Kernel::translationComponents]; for(int c=0;c<Kernel::translationComponents;++c) { data[c]=cache((id*P+j)*P+a,c); }
                                    kernel.applyTranslation(data,q,sum);
                                }
                                else { kernel.add(r,q,sum); }
                            }
                            else { kernel.add(r,q,sum); }
                        }
                    }
                    for(int c=0;c<FieldComponents;++c) { L(b*P+j,c)+=sum[c]; }
                });
                if(l==depth_) { continue; }
                const auto childStart=levelStart_[l+1]; const int m=2*n;
                Kokkos::parallel_for("FMM L2L",Kokkos::RangePolicy<ExecutionSpace>(ex,0,m*m*m*P),KOKKOS_LAMBDA(const int k)
                {
                    const int b=k/P,j=k%P,x=b%m,y=(b/m)%m,z=b/(m*m); if(!tc(childStart+b)) { return; }
                    const auto parent=start+x/2+n*(y/2+n*(z/2)); const int oct=(x&1)+2*(y&1)+4*(z&1);
                    double sum[FieldComponents]={};
                    for(int a=0;a<P;++a) { const double w=tr(oct*P+a,j); for(int c=0;c<FieldComponents;++c) { sum[c]+=w*L(parent*P+a,c); } }
                    for(int c=0;c<FieldComponents;++c) { L((childStart+b)*P+j,c)=sum[c]; }
                });
            }
            if(interpolateNearBoxes)
            {
                // Smooth-kernel self and neighbor translations. Each target leaf
                // owns its locals, so no atomics or particle-pair work is needed.
                Kokkos::parallel_for("FMM smooth near M2L",Kokkos::RangePolicy<ExecutionSpace>(ex,0,leaves*P),KOKKOS_LAMBDA(const int k)
                {
                    const auto b=k/P;const int j=k%P;
                    if(!tc(leaf+b)) {
              return;
            }
            double sum[FieldComponents]={};
                    for(std::size_t pair=no(b);pair<no(b+1);++pair)
                    {
                        const auto sb=ni(pair);
                        for(int a=0;a<P;++a)
                        {
                            double r[3],q[SourceComponents];
                            for(int c=0;c<3;++c) {
                  r[c]=xyz((leaf+b)*P+j,c)-xyz((leaf+sb)*P+a,c);
                }
                for(int c=0;c<SourceComponents;++c) {
                  q[c]=M((leaf+sb)*P+a,c);
                }
                kernel.add(r,q,sum);
                        }
                    }
                    for(int c=0;c<FieldComponents;++c) {
              L((leaf+b)*P+j,c)+=sum[c];
            }
          });
            }
            checkpoint(2);
            Kokkos::parallel_for("FMM L2P and near",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex,0,nt),KOKKOS_LAMBDA(const std::size_t t)
            {
                const auto b=tl(t); double sum[FieldComponents]={};
                for(int j=0;j<P;++j)
                {
                    const double w=tw(t,j%p)*tw(t,p+(j/p)%p)*tw(t,2*p+j/(p*p));
                    for(int c=0;c<FieldComponents;++c) { sum[c]+=w*L((leaf+b)*P+j,c); }
                }
                for(std::size_t pair=no(b);!interpolateNearBoxes && !externalNear && pair<no(b+1);++pair)
                {
                    const auto sb=ni(pair);
                    if(interpolateNear && so(sb+1)-so(sb)>std::size_t(P))
                    {
                        // Source interpolation only: retain exact target locations and avoid
                        // another interpolation across the near-field variation.
                        for(int a=0;a<P;++a)
                        {
                            double r[3],q[SourceComponents];
                            for(int c=0;c<3;++c) { r[c]=tgt(t,c)-xyz((leaf+sb)*P+a,c); }
                            for(int c=0;c<SourceComponents;++c) { q[c]=M((leaf+sb)*P+a,c); }
                            kernel.add(r,q,sum);
                        }
                        continue;
                    }
                    for(std::size_t a=so(sb);a<so(sb+1);++a)
                    {
                        const auto s=si(a); double r[3],q[SourceComponents]; for(int c=0;c<3;++c) { r[c]=tgt(t,c)-src(s,c); }
                        for(int c=0;c<SourceComponents;++c) { q[c]=tensor(s,c); } kernel.add(r,q,sum);
                    }
                }
                OutputPolicy::store(output,t,sum);
            });
            checkpoint(3);
            // Work arrays are reused; callers may consume the result after this fence.
            ex.fence();
        }
    };
}
#endif
