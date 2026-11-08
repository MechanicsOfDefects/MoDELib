/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_KokkosParticleFMM_H_
#define model_KokkosParticleFMM_H_
#include <cfloat>
#include <Eigen/Dense>
#include <cstdlib>
#include <Kokkos_Core.hpp>
#include <FMM/BoxParticleIndex.h>
#include <FMM/ShiftedBoxNeighbors.h>
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
        View imageShifts_;
        Indices farImages_,nearImages_;
        std::vector<std::size_t> hostTranslationImages_;
        View centers_,nodeCoordinates_,sourceWeights_,targetWeights_,transfer_,moments_,locals_;
        Indices sourceOffsets_,sourceIDs_,targetLeaves_,m2lOffsets_,m2lSources_,nearOffsets_,nearSources_,sourceCounts_,targetCounts_;
        std::size_t farPairs_=0,nearPairs_=0;
        bool externalNear_=false,interpolateNear_=false,interpolateNearBoxes_=false,profile_=false;
        std::array<double,3> translationTimings_{}; // M2L, L2L, smooth near M2L (ms)
        std::array<double,4> lastTimings_{}; // cache preparation, upward, downward, target/near (ms)
        std::size_t cacheBudget_=256u*1024u*1024u;
        View translations_;
        View integratedLeafMoments_;
        Indices translationIDs_,nearTranslationIDs_;
        std::vector<std::size_t> hostTranslationIDs_,hostNearTranslationIDs_,hostTranslationTargets_,hostTranslationSources_;
        std::size_t farTranslationClasses_=0;
        std::vector<double> translationKey_{};
        std::type_index translationType_{typeid(void)};
        bool cacheReady_=false;
        // Diagnostic reference switch; ordinary runs use grouped M2L.
        bool groupedM2L_=std::getenv("MODELIB_FMM_SCALAR_M2L")==nullptr;
        bool scalarGroupedMatrix_=false;
        // Immutable geometry grouping; numeric packing is refreshed for each RHS.
        std::vector<std::size_t> hostPairSources_;
        std::vector<std::size_t> groupLevelOffsets_;
        Indices groupOffsets_,groupCacheIDs_,groupMatrixIDs_,groupSources_,pairGroupSlots_;
        View groupInputs_,groupOutputs_,groupMatrices_;
        std::size_t groupedPairs_=0;
        template<class Kernel> void prepareGroups(const ExecutionSpace& ex,const Kernel& kernel,
                                                  const std::vector<std::size_t>& cacheIDs)
        {
            const auto missing=std::numeric_limits<std::size_t>::max();
            std::vector<std::size_t> slots(farPairs_,missing),offsets{0},ids,matrixIDs,sources;
            groupLevelOffsets_.assign(depth_+2,0);
            groupedPairs_=0;
            groupMatrices_=View();
            scalarGroupedMatrix_=SourceComponents==1 && FieldComponents==1 && Kernel::translationComponents==1;
            // Eigen GEMM is host-only. Device execution retains the original kernel.
            if constexpr(std::is_same_v<typename ExecutionSpace::memory_space,Kokkos::HostSpace>)
            {
                if(groupedM2L_)
                {
                    std::vector<std::vector<std::size_t>> groups(hostTranslationTargets_.size());
                    for(std::size_t pair=0;pair<cacheIDs.size();++pair)
                    {
                        if(cacheIDs[pair]!=missing)
                        {
                            groups[hostTranslationIDs_[pair]].push_back(pair);
                        }
                    }
                    const std::size_t matrixEntries=std::size_t(nodes_)*nodes_*SourceComponents*FieldComponents;
                    const auto used=translations_.size()*sizeof(double);
                    const auto matrixBudget=cacheBudget_>used?(cacheBudget_-used)/(matrixEntries*sizeof(double)):0;
                    std::size_t expanded=0;
                    for(unsigned level=0;level<=depth_;++level)
                    {
                        groupLevelOffsets_[level]=ids.size();
                        const auto end=levelStart_[level]+(std::size_t(1)<<(3*level));
                        for(std::size_t id=0;id<groups.size();++id)
                        {
                            const auto& pairs=groups[id];
                            if(pairs.size()<4 || hostTranslationTargets_[id]<levelStart_[level] || hostTranslationTargets_[id]>=end)
                            {
                                continue;
                            }
                            if constexpr(SourceComponents!=1 || FieldComponents!=1 || Kernel::translationComponents!=1)
                            {
                                if(expanded>=matrixBudget)
                                {
                                    continue;
                                }
                            }
                            ids.push_back(cacheIDs[pairs.front()]);
                            matrixIDs.push_back(expanded);
                            if constexpr(SourceComponents!=1 || FieldComponents!=1 || Kernel::translationComponents!=1)
                            {
                                ++expanded;
                            }
                            for(const auto pair:pairs)
                            {
                                slots[pair]=sources.size();
                                sources.push_back(hostPairSources_[pair]);
                            }
                            offsets.push_back(sources.size());
                        }
                    }
                    groupLevelOffsets_[depth_+1]=ids.size();
                    groupMatrices_=View("FMM expanded grouped translations",expanded*matrixEntries,1);
                    if constexpr(SourceComponents!=1 || FieldComponents!=1 || Kernel::translationComponents!=1)
                    {
                        const auto cache=translations_,matrices=groupMatrices_;
                        const auto cacheIndices=upload(ex,"FMM expanded translation IDs",ids);
                        const int P=nodes_,rows=P*FieldComponents,cols=P*SourceComponents;
                        Kokkos::parallel_for("FMM expand grouped translations",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex,0,expanded*matrixEntries),KOKKOS_LAMBDA(const std::size_t k)
                        {
                            const auto group=k/matrixEntries,entry=k%matrixEntries;
                            const int row=entry%rows,col=entry/rows,j=row/FieldComponents,a=col/SourceComponents;
                            double data[Kernel::translationComponents],q[SourceComponents]={},out[FieldComponents]={};
                            for(int c=0;c<Kernel::translationComponents;++c)
                            {
                                data[c]=cache((cacheIndices(group)*P+j)*P+a,c);
                            }
                            q[col%SourceComponents]=1.0;
                            kernel.applyTranslation(data,q,out);
                            matrices(k,0)=out[row%FieldComponents];
                        });
                    }
                }
            }
            groupedPairs_=sources.size();
            groupOffsets_=upload(ex,"FMM group offsets",offsets);
            groupCacheIDs_=upload(ex,"FMM group cache IDs",ids);
            groupMatrixIDs_=upload(ex,"FMM group matrix IDs",matrixIDs);
            groupSources_=upload(ex,"FMM group sources",sources);
            pairGroupSlots_=upload(ex,"FMM pair group slots",slots);
            groupInputs_=View("FMM grouped inputs",nodes_*SourceComponents,sources.size());
            groupOutputs_=View("FMM grouped outputs",nodes_*FieldComponents,sources.size());
            ex.fence();
        }
        void applyGroups(const ExecutionSpace& ex,unsigned level)
        {
            if constexpr(std::is_same_v<typename ExecutionSpace::memory_space,Kokkos::HostSpace>)
            {
                if(groupLevelOffsets_.empty())
                {
                    return;
                }
                const auto first=groupLevelOffsets_[level],last=groupLevelOffsets_[level+1];
                const auto offsets=groupOffsets_,sources=groupSources_,cacheIDs=groupCacheIDs_,matrixIDs=groupMatrixIDs_;
                const auto input=groupInputs_,output=groupOutputs_,moments=moments_,cache=translations_,matrices=groupMatrices_;
                const int P=nodes_,rows=P*FieldComponents,cols=P*SourceComponents;
                const bool scalar=scalarGroupedMatrix_;
                // Each task owns its packed columns and output. Scatter is performed
                // by the target-owned M2L loop, preserving pair order without atomics.
                Kokkos::parallel_for("FMM grouped GEMM",Kokkos::RangePolicy<ExecutionSpace>(ex,first,last),[=](const int group)
                {
                    const auto begin=offsets(group),end=offsets(group+1);
                    for(auto i=begin;i<end;++i)
                    {
                        for(int a=0;a<P;++a)
                        {
                            for(int c=0;c<SourceComponents;++c)
                            {
                                input(a*SourceComponents+c,i)=moments(sources(i)*P+a,c);
                            }
                        }
                    }
                    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
                    const Eigen::Map<const Matrix> Q(input.data()+begin*cols,cols,end-begin);
                    Eigen::Map<Matrix> Y(output.data()+begin*rows,rows,end-begin);
                    if(scalar)
                    {
                        using RowMatrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
                        const Eigen::Map<const RowMatrix> T(cache.data()+cacheIDs(group)*P*P,P,P);
                        Y.noalias()=T*Q;
                    }
                    else
                    {
                        const Eigen::Map<const Matrix> T(matrices.data()+matrixIDs(group)*std::size_t(rows)*cols,rows,cols);
                        Y.noalias()=T*Q;
                    }
                });
            }
        }
        template<class Kernel> void prepareTranslations(const ExecutionSpace& ex,const Kernel& kernel)
        {
            if constexpr(HasFMMTranslationCache<Kernel>::value)
            {
                const auto rawKey=kernel.translationKey();
                const std::vector<double> key(rawKey.begin(),rawKey.end());
                if(cacheReady_ && translationType_==std::type_index(typeid(Kernel)) && translationKey_==key) { return; }
                ex.fence();
                const std::size_t perClass=std::size_t(nodes_)*nodes_;
                // Prioritize the repeatedly used self/neighbor operators within the
                // same memory budget. Unsplit plans allocate only far operators.
                const auto total=interpolateNearBoxes_?hostTranslationTargets_.size():farTranslationClasses_;
                const auto count=std::min(total,cacheBudget_/(perClass*Kernel::translationComponents*sizeof(double)));
                const auto missing=std::numeric_limits<std::size_t>::max();
                std::vector<std::size_t> remap(hostTranslationTargets_.size(),missing),selectedTargets,selectedSources,selectedImages;
                const auto select=[&](std::size_t begin,std::size_t end)
                {
                    for(auto id=begin;id<end && selectedTargets.size()<count;++id)
                    {
                        remap[id]=selectedTargets.size();
                        selectedTargets.push_back(hostTranslationTargets_[id]);
                        selectedSources.push_back(hostTranslationSources_[id]);
                        selectedImages.push_back(hostTranslationImages_[id]);
                    }
                };
                if(interpolateNearBoxes_)
                {
                    select(farTranslationClasses_,hostTranslationTargets_.size());
                }
                select(0,farTranslationClasses_);
                auto ids=hostTranslationIDs_,nearIDs=hostNearTranslationIDs_;
                for(auto& id:ids) { id=remap[id]; }
                for(auto& id:nearIDs) { id=remap[id]; }
                translationIDs_=upload(ex,"FMM translation IDs",ids);
                nearTranslationIDs_=upload(ex,"FMM near translation IDs",nearIDs);
                translations_=View("FMM cached translation coefficients",count*perClass,Kernel::translationComponents);
                const auto xyz=nodeCoordinates_,cache=translations_,shifts=imageShifts_;
                const auto images=upload(ex,"FMM translation images",selectedImages);
                const auto targets=upload(ex,"FMM translation targets",selectedTargets);
                const auto sources=upload(ex,"FMM translation sources",selectedSources);
                const int P=nodes_;
                Kokkos::parallel_for("FMM prepare translations",Kokkos::RangePolicy<ExecutionSpace,Kokkos::IndexType<std::size_t>>(ex,0,count*perClass),KOKKOS_LAMBDA(const std::size_t k)
                {
                    const auto id=k/(P*P),j=(k/P)%P,a=k%P;
                    double r[3],data[Kernel::translationComponents];
                    for(int c=0;c<3;++c) { r[c]=xyz(targets(id)*P+j,c)-xyz(sources(id)*P+a,c)-shifts(images(id),c); }
                    kernel.prepareTranslation(r,data);
                    for(int c=0;c<Kernel::translationComponents;++c) { cache(k,c)=data[c]; }
                });
                prepareGroups(ex,kernel,ids);
                ex.fence(); translationKey_=key; translationType_=std::type_index(typeid(Kernel)); cacheReady_=true;
            }
        }
        template<class Kernel> KOKKOS_INLINE_FUNCTION static void applyTranslation(
            const Kernel& kernel,const View& cache,std::size_t id,int P,int j,int a,
            const double* r,const double* q,double* sum)
        {
            if constexpr(HasFMMTranslationCache<Kernel>::value)
            {
                if(id!=std::numeric_limits<std::size_t>::max())
                {
                    double data[Kernel::translationComponents];
                    for(int c=0;c<Kernel::translationComponents;++c)
                    {
                        data[c]=cache((id*P+j)*P+a,c);
                    }
                    kernel.applyTranslation(data,q,sum);
                    return;
                }
            }
            kernel.add(r,q,sum);
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
        KokkosParticleFMM(const ExecutionSpace& ex,const BoxTree<3>& tree,const Positions& sources,const Positions& targets,int order=4,const std::vector<Eigen::Vector3d>& shifts={Eigen::Vector3d::Zero()})
            :tree_(&tree),revision_(tree.revision()),depth_(tree.depth()),order_(checkedOrder(order)),nodes_(order_*order_*order_),ns_(sources[0].size()),nt_(targets[0].size())
        {
            validateImageShifts(shifts);
            imageShifts_=View("FMM image shifts",shifts.size(),3);
            auto imageHost=Kokkos::create_mirror(imageShifts_);
            for(std::size_t i=0;i<shifts.size();++i) { for(int c=0;c<3;++c) { imageHost(i,c)=shifts[i][c]; } }
            Kokkos::deep_copy(ex,imageShifts_,imageHost);
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
                const int n=1<<l; const Eigen::Matrix3d cellBasis=tree.basis()/n; diameter[l]=tree.cellDiameter(l);
                for(int z=0;z<n;++z) { for(int y=0;y<n;++y) { for(int x=0;x<n;++x)
                {
                    const auto b=levelStart_[l]+x+n*(y+n*z); center[b]=tree.physical(Eigen::Vector3d(x+.5,y+.5,z+.5)/n);
                    for(int a=0;a<3;++a) { ch(b,a)=center[b][a]; }
                    for(int j=0;j<nodes_;++j)
                    {
                        const int ij[3]={j%order,(j/order)%order,j/(order*order)};
                        const Eigen::Vector3d position=center[b]+.5*cellBasis*Eigen::Vector3d(nodes[ij[0]],nodes[ij[1]],nodes[ij[2]]);
                        for(int a=0;a<3;++a) { nh(b*nodes_+j,a)=position[a]; }
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
            const Eigen::Matrix3d leafInverse=double(1u<<depth_)*tree.inverseBasis();
            const auto setWeights=[&](const Positions& pts,View view,const std::vector<std::size_t>& leafOfPoint)
            {
                auto host=Kokkos::create_mirror(view);
                for(std::size_t p=0;p<pts[0].size();++p) { for(int a=0;a<3;++a)
                {
                    const Eigen::Vector3d q=2*leafInverse*(Eigen::Vector3d(pts[0][p],pts[1][p],pts[2][p])-center[leafStart_+leafOfPoint[p]]);
                    const auto w=weights(q[a],nodes);
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
            std::vector<std::vector<std::size_t>> far(boxCount_),near(leaves),farImage(boxCount_),nearImage(leaves);
            std::size_t currentImage=0;
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
                if((center[t]-center[s]-shifts[currentImage]).norm()>1.01*diameter[l]) { far[t].push_back(s); farImage[t].push_back(currentImage); ++farPairs_; return; }
                if(l==depth_) { near[t-leafStart_].push_back(s-leafStart_); nearImage[t-leafStart_].push_back(currentImage); ++nearPairs_; return; }
                for(int a=0;a<8;++a) { for(int b=0;b<8;++b) { walk(l+1,childID(l,t,a),childID(l,s,b)); } }
            };
            for(currentImage=0;currentImage<shifts.size();++currentImage) { walk(0,0,0); }
            const auto pack=[&](const auto& rows,const char* ol,const char* il,Indices& offsets,Indices& ids)
            {
                std::vector<std::size_t> o(rows.size()+1),v;
                for(std::size_t b=0;b<rows.size();++b) { v.insert(v.end(),rows[b].begin(),rows[b].end()); o[b+1]=v.size(); }
                offsets=upload(ex,ol,o); ids=upload(ex,il,v);
            };
            pack(far,"FMM M2L offsets","FMM M2L sources",m2lOffsets_,m2lSources_);
            pack(near,"FMM near offsets","FMM near sources",nearOffsets_,nearSources_);
            std::vector<std::size_t> fiHost,niHost;
            for(const auto& row:farImage) { fiHost.insert(fiHost.end(),row.begin(),row.end()); }
            for(const auto& row:nearImage) { niHost.insert(niHost.end(),row.begin(),row.end()); }
            farImages_=upload(ex,"FMM far images",fiHost);nearImages_=upload(ex,"FMM near images",niHost);
            // Translation invariance: share one operator for each level and integer displacement.
            std::map<std::array<int,5>,std::size_t> classes;
            for(unsigned l=0;l<=depth_;++l)
            {
                const int n=1<<l;
                for(std::size_t b=levelStart_[l];b<levelStart_[l]+std::size_t(n)*n*n;++b)
                {
                    const int t=static_cast<int>(b-levelStart_[l]);
                    for(std::size_t pair=0;pair<far[b].size();++pair)
                    {
                        const auto source=far[b][pair],image=farImage[b][pair];
                        const int s=static_cast<int>(source-levelStart_[l]);
                        const std::array<int,5> key{static_cast<int>(l),t%n-s%n,(t/n)%n-(s/n)%n,t/(n*n)-s/(n*n),static_cast<int>(image)};
                        const auto result=classes.emplace(key,classes.size());
                        if(result.second) { hostTranslationTargets_.push_back(b); hostTranslationSources_.push_back(source); hostTranslationImages_.push_back(image); }
                        hostTranslationIDs_.push_back(result.first->second);
                        hostPairSources_.push_back(source);
                    }
                }
            }
            farTranslationClasses_=classes.size();
            const int leafN=1<<depth_;
            for(std::size_t b=0;b<near.size();++b)
            {
                const int t=static_cast<int>(b);
                for(std::size_t pair=0;pair<near[b].size();++pair)
                {
                    const auto source=near[b][pair],image=nearImage[b][pair];
                    const int si=static_cast<int>(source);
                    const std::array<int,5> key{static_cast<int>(depth_),t%leafN-si%leafN,(t/leafN)%leafN-(si/leafN)%leafN,t/(leafN*leafN)-si/(leafN*leafN),static_cast<int>(image)};
                    const auto result=classes.emplace(key,classes.size());
                    if(result.second)
                    {
                        hostTranslationTargets_.push_back(leafStart_+b);
                        hostTranslationSources_.push_back(leafStart_+source);
                        hostTranslationImages_.push_back(image);
                    }
                    hostNearTranslationIDs_.push_back(result.first->second);
                }
            }
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
        void setGroupedM2L(bool enabled)
        {
            if(groupedM2L_!=enabled)
            {
                groupedM2L_=enabled;
                cacheReady_=false;
            }
        }
        std::size_t groupedPairCount() const { return groupedPairs_; }
        std::size_t groupedMatrixBytes() const { return groupMatrices_.size()*sizeof(double); }
        void setProfiling(bool enabled) { profile_=enabled; }
        const std::array<double,4>& lastTimings() const { return lastTimings_; }
        const std::array<double,3>& translationTimings() const { return translationTimings_; }
        void setExternalNear(bool enabled){externalNear_=enabled;}
        void setInterpolateNearBoxes(bool enabled) { if(interpolateNearBoxes_!=enabled) { interpolateNearBoxes_=enabled; cacheReady_=false; } }
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
            if constexpr(!HasFMMTranslationCache<Kernel>::value)
            {
                pairGroupSlots_=Indices();
                groupLevelOffsets_.clear();
                cacheReady_=false;
            }
            translationTimings_.fill(0.0);
            auto translationStamp=stamp;
            const auto translationCheckpoint=[&](int stage)
            {
                if(profile_)
                {
                    ex.fence();
                    const auto now=Clock::now();
                    translationTimings_[stage]+=std::chrono::duration<double,std::milli>(now-translationStamp).count();
                    translationStamp=now;
                }
            };
            prepareTranslations(ex,kernel); checkpoint(0);
            const bool externalNear=externalNear_,interpolateNear=interpolateNear_,interpolateNearBoxes=interpolateNearBoxes_;
            const auto cache=translations_; const auto cacheIDs=translationIDs_,nearCacheIDs=nearTranslationIDs_;
            const int p=order_,P=nodes_; const auto leaf=leafStart_;
            const auto M=moments_,L=locals_,sw=sourceWeights_,tw=targetWeights_,tr=transfer_,xyz=nodeCoordinates_;
            const auto so=sourceOffsets_,si=sourceIDs_,tl=targetLeaves_,sc=sourceCounts_,tc=targetCounts_;
            const auto groupSlots=pairGroupSlots_; const auto groupOutput=groupOutputs_;
            const auto fo=m2lOffsets_,fi=m2lSources_,no=nearOffsets_,ni=nearSources_;
            const auto shifts=imageShifts_;const auto fimages=farImages_,nimages=nearImages_;
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
            translationStamp=Clock::now();
            // Shifted source roots can themselves be admissible M2L pairs.
            for(unsigned l=0;l<=depth_;++l)
            {
                const auto start=levelStart_[l]; const int n=1<<l;
                applyGroups(ex,l);
                Kokkos::parallel_for("FMM M2L",Kokkos::RangePolicy<ExecutionSpace>(ex,0,n*n*n*P),KOKKOS_LAMBDA(const int k)
                {
                    const auto b=start+k/P; const int j=k%P; if(!tc(b)) { return; }
                    double sum[FieldComponents]={};
                    for(std::size_t pair=fo(b);pair<fo(b+1);++pair)
                    {
                        if(groupSlots.extent(0) && groupSlots(pair)!=std::numeric_limits<std::size_t>::max())
                        {
                            for(int c=0;c<FieldComponents;++c)
                            {
                                sum[c]+=groupOutput(j*FieldComponents+c,groupSlots(pair));
                            }
                            continue;
                        }
                        const auto s=fi(pair);
                        for(int a=0;a<P;++a)
                        {
                            double r[3],q[SourceComponents]; for(int c=0;c<3;++c) { r[c]=xyz(b*P+j,c)-xyz(s*P+a,c)-shifts(fimages(pair),c); }
                            for(int c=0;c<SourceComponents;++c) { q[c]=M(s*P+a,c); }
                            std::size_t id=std::numeric_limits<std::size_t>::max();
                            if constexpr(HasFMMTranslationCache<Kernel>::value) { id=cacheIDs(pair); }
                            applyTranslation(kernel,cache,id,P,j,a,r,q,sum);
                        }
                    }
                    for(int c=0;c<FieldComponents;++c) { L(b*P+j,c)+=sum[c]; }
                });
                translationCheckpoint(0);
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
                translationCheckpoint(1);
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
                  r[c]=xyz((leaf+b)*P+j,c)-xyz((leaf+sb)*P+a,c)-shifts(nimages(pair),c);
                }
                for(int c=0;c<SourceComponents;++c) {
                  q[c]=M((leaf+sb)*P+a,c);
                }
                std::size_t id=std::numeric_limits<std::size_t>::max();
                if constexpr(HasFMMTranslationCache<Kernel>::value) { id=nearCacheIDs(pair); }
                applyTranslation(kernel,cache,id,P,j,a,r,q,sum);
                        }
                    }
                    for(int c=0;c<FieldComponents;++c) {
              L((leaf+b)*P+j,c)+=sum[c];
            }
          });
            }
            translationCheckpoint(2);
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
                            for(int c=0;c<3;++c) { r[c]=tgt(t,c)-xyz((leaf+sb)*P+a,c)-shifts(nimages(pair),c); }
                            for(int c=0;c<SourceComponents;++c) { q[c]=M((leaf+sb)*P+a,c); }
                            kernel.add(r,q,sum);
                        }
                        continue;
                    }
                    for(std::size_t a=so(sb);a<so(sb+1);++a)
                    {
                        const auto s=si(a); double r[3],q[SourceComponents]; for(int c=0;c<3;++c) { r[c]=tgt(t,c)-src(s,c)-shifts(nimages(pair),c); }
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
