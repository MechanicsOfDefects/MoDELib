#include <cfloat>
#include <FMM/BoxBoundaryTolerance.h>
#include <Utilities/Profiling.h>
#include <DislocationDynamics/StraightDislocationSegmentFactory.h>
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
#include <DislocationDynamics/StraightDislocationSegmentCai.h>
#include <iostream>
#include <DislocationDynamics/StraightDislocationSegmentHelmholtz.h>
#include <FMM/KokkosFMMInteraction.h>
#include <FMM/HelmholtzDiffusionKernel.h>
#include <DislocationDynamics/SegmentConcentrationPlan.h>
#include <Utilities/SolverStageTimer.h>
#include <FMM/HelmholtzStressKernelFamily.h>
#include <FMM/CaiStressKernelFamily.h>
#endif
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <atomic>
#include <Utilities/ParallelFor.h>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <limits>
namespace model
{
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
    namespace
    {
        using Ex=Kokkos::DefaultExecutionSpace;
        using View=Kokkos::View<double**,Kokkos::LayoutLeft,Ex>;
        using Indices=Kokkos::View<std::size_t**,Kokkos::LayoutLeft,Ex>;
        using Positions=std::array<std::vector<double>,3>;
        struct Registrations
        {
            std::shared_ptr<FMMManager> manager;
            std::string prefix;
            std::vector<std::string> kernels,interactions;
            explicit Registrations(std::shared_ptr<FMMManager> m):manager(std::move(m))
            {
                if(!manager) {
      throw std::invalid_argument("Missing shared FMM manager");
    }
    if(!Kokkos::is_initialized()||Kokkos::is_finalized()) {
      throw std::logic_error("Initialize Kokkos before screened assembly");
    }
    static std::atomic<std::size_t> serial{0};
                prefix="DD/screened/"+std::to_string(serial.fetch_add(1))+"/";
            }
            ~Registrations()
            {
                for(const auto& n:interactions) {
      manager->unregisterInteraction(n);
    }
    for(const auto& n:kernels) {
      manager->unregisterKernel(n);
    }
  }
            template<class Family> auto addKernel(const std::string& name,Family family)
            {
                const auto full=prefix+name; kernels.push_back(full);
                return manager->registerKernel(full,std::move(family));
            }
            template<class Interaction,class Handle,class... Args> auto addInteraction(const std::string& name,const Handle& h,Args&&... args)
            {
                const auto full=prefix+name; interactions.push_back(full);
                return manager->registerInteraction<Interaction>(full,h,std::forward<Args>(args)...);
            }
        };
        FMMInteractionOptions interactionOptions(const ScreenedDislocationOptions& options)
        {
            FMMInteractionOptions opt;
            opt.localImageShifts=options.localImageShifts;
            opt.boxShortRange=options.boxShortRange;opt.interpolateScreenedNear=options.interpolateScreenedNear;opt.interpolationOrder=options.interpolationOrder;opt.shortCutoff=options.shortCutoff;
            opt.translationCacheBytes=options.translationCacheBytes;return opt;
        }
        Positions positions(const std::vector<ClimbQuadratureMap>& points,std::size_t nodes)
        {
            Positions p;
            for(const auto& q:points)
            {
                for(int j=0;j<2;++j) {
      if(q.nodes[j]>=nodes||!std::isfinite(q.weights[j])) {
        throw std::invalid_argument("Invalid climb quadrature map");
      }
    }
    for(int c=0;c<3;++c) { if(!std::isfinite(q.position[c])) {
        throw std::invalid_argument("Invalid climb quadrature position");
      }
      p[c].push_back(q.position[c]); }
            }
            return p;
        }
        void uploadMap(const Ex& ex,const std::vector<ClimbQuadratureMap>& points,Indices& ids,View& weights)
        {
            ids=Indices("climb endpoint indices",points.size(),2);weights=View("climb endpoint weights",points.size(),2);
            auto ih=Kokkos::create_mirror(ids);auto wh=Kokkos::create_mirror(weights);
            for(std::size_t q=0;q<points.size();++q) {
    for(int j=0;j<2;++j){ih(q,j)=points[q].nodes[j];wh(q,j)=points[q].weights[j];}
  }
  Kokkos::deep_copy(ex,ids,ih);Kokkos::deep_copy(ex,weights,wh);ex.fence();
        }
    }
    struct FMMStressEvaluator::Impl
    {
        virtual ~Impl()=default;
        virtual DislocationStressArrays<3> evaluate()=0;
        virtual Eigen::Matrix3d stressAt(const Eigen::Vector3d&) const=0;
    };
    template<class Family> struct FMMStressEvaluator::StressImpl : Impl
    {
        Registrations registrations;
        Ex ex;
        double mu_,nu_;
        using Interaction=KokkosFMMInteraction<Family>;
        std::shared_ptr<Interaction> interaction;
        View strengths,output;
        std::shared_ptr<const BoxTree<3>> tree_;
        std::uint64_t revision_;
        DislocationSourceArrays<3> sources_;
        BoxParticleIndex<3> index_;
        Family family_;
        ScreenedDislocationOptions options_;
        std::size_t originalTargets_;
        std::vector<double> leafLocals_;
        std::vector<std::size_t> lexToBox_;
        std::shared_ptr<const StraightDislocationSegmentBase> segmentKernel_;
        std::vector<std::vector<helmholtzSegment::Segment>> segmentBins_;
        std::unique_ptr<ShiftedBoxNeighbors> nearNeighborhood_;
        DislocationTargetArrays<3> targetSnapshot_;
        StressImpl(std::shared_ptr<FMMManager> m,std::shared_ptr<const BoxTree<3>> tree,
             const DislocationSourceArrays<3>& src,const DislocationTargetArrays<3>& tgt,
             double mu,double nu,const ScreenedDislocationOptions& opt,const std::vector<helmholtzSegment::Segment>& segments):registrations(std::move(m)),
             mu_(mu),nu_(nu),strengths("stress quadrature strengths",src.size(),9),tree_(tree),revision_(tree?tree->revision():0),
             sources_(src),family_(mu,nu,opt.coreWidth,opt.screenWidth,opt.screened),options_(opt),originalTargets_(tgt.size())
        {
            if(!tree_) {
      throw std::invalid_argument("Missing stress tree");
    }
    targetSnapshot_=tgt;
            nearNeighborhood_=std::make_unique<ShiftedBoxNeighbors>(*tree_,options_.localImageShifts);
            if(opt.segmentShortRange){
                if(!opt.boxShortRange||(opt.screened&&(!opt.interpolateScreenedNear||(opt.coreType==CoreType::Cai)))) {
        throw std::invalid_argument("Segment short stress requires box-local screened Helmholtz FMM");
      }
      SegmentCoreParameters parameters;
      parameters.coreWidth=opt.coreWidth;parameters.ewaldWidth=opt.screenWidth;
      parameters.mu=mu;parameters.nu=nu;parameters.order=opt.segmentChebyshevOrder;
      parameters.resolution=opt.segmentTableResolution;parameters.maxDistance=2*tree_->cellDiameter(tree_->depth());
      segmentKernel_=makeStraightDislocationSegmentCore(opt.coreType,opt.screened,parameters);
                if(!options_.segmentSplits)
                {
                    std::vector<SegmentSplitSnapshot::Endpoints> geometry;
                    for(const auto& segment:segments)
                    {
                        geometry.emplace_back(segment.p0,segment.p1);
                    }
                    options_.segmentSplits=std::make_shared<SegmentSplitSnapshot>(*tree_,SegmentSplitSnapshot::geometry(geometry));
                }
                if(!options_.segmentSplits->matches(*tree_))
                {
                    throw std::logic_error("Stale shared stress segment split tree");
                }
                segmentBins_.resize(options_.segmentSplits->leafCount());
                for(const auto& segment:segments)
                {
                    for(const auto& piece:options_.segmentSplits->pieces(segment.p0,segment.p1))
                    {
                        segmentBins_[piece.leaf].push_back({piece.p0,piece.p1,segment.b});
                    }
                }
                std::size_t count=0;for(const auto& bin:segmentBins_) {
        count+=bin.size();
      }
      std::cout<<(opt.screened?"  Ewald short stress: ":"  FMM full segment near stress: ")<<segments.size()<<" segments, "<<count<<" box pieces; Chebyshev order="<<opt.segmentChebyshevOrder<<", resolution="<<opt.segmentTableResolution<<std::endl;
            }
            if(options_.segmentShortRange){
                // Mark occupied leaves only; source coefficients are integrated below.
                sources_=DislocationSourceArrays<3>();
                for(std::size_t b=0;b<segmentBins_.size();++b) {
        if(!segmentBins_[b].empty()){
                    const auto point=tree_->cellCenter(tree_->depth(),b);
                    std::array<double,3> x{point[0],point[1],point[2]};
                    sources_.positions.append(x);sources_.weightedBurgersTensors.append(std::array<double,9>{});
                }
      }
      strengths=View("stress occupied leaf markers",sources_.size(),9);
            }
            const auto& effective=sources_;
            auto targets=tgt.positions.values;
            // Occupy every leaf so cached locals support arbitrary FEM and DD targets.
            const int n=1<<tree_->depth();const Eigen::Vector3d width=Eigen::Vector3d::Constant(1.0/n);
            for(int z=0;z<n;++z) {
      for(int y=0;y<n;++y) {
        for(int x=0;x<n;++x)
            {const Eigen::Vector3d point=tree_->cellCenter(tree_->depth(),x+n*(y+n*z));for(int c=0;c<3;++c) {
            targets[c].push_back(point[c]);
          }
        }
      }
    }
    index_.rebuild(*tree_,effective.positions.values,targets);
            lexToBox_.resize(n*n*n);
            for(std::size_t i=0;i<lexToBox_.size();++i){const auto key=tree_->boxes(tree_->depth())[i]->key();const auto w=BoxTree<3>::latticeSize/n;lexToBox_[key[0]/w+n*(key[1]/w+n*(key[2]/w))]=i;}
            output=View("stress quadrature outputs",targets[0].size(),6);
            auto kernel=registrations.addKernel("stress",family_);
            auto io=interactionOptions(opt);if(options_.segmentShortRange){io.buildShortNeighbors=false;io.externalNear=!opt.screened;}
            interaction=registrations.addInteraction<Interaction>("stress/interaction",kernel,ex,std::move(tree),effective.positions.values,targets,io);
            auto host=Kokkos::create_mirror(strengths);
            for(int c=0;c<9;++c)
            {
                if(effective.weightedBurgersTensors.values[c].size()!=effective.size()) {
        throw std::invalid_argument("Invalid stress source strengths");
      }
      for(std::size_t q=0;q<effective.size();++q){const auto value=effective.weightedBurgersTensors.values[c][q];if(!std::isfinite(value)) {
          throw std::invalid_argument("Nonfinite stress source");
        }
        host(q,c)=value;}
            }
            Kokkos::deep_copy(ex,strengths,host);ex.fence();
            if(options_.segmentShortRange){
                const int p=opt.interpolationOrder,P=p*p*p;
                View coefficients("segment integrated stress moments",segmentBins_.size()*P,9);
                auto h=Kokkos::create_mirror(coefficients);Kokkos::deep_copy(h,0.0);
                const SegmentMomentIntegrator helper(p,tree_->inverseBasis());
#ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for(std::size_t b=0;b<segmentBins_.size();++b){
                    const Eigen::Vector3d center=tree_->cellCenter(tree_->depth(),b);
                    for(const auto& piece:segmentBins_[b]){
                        const auto values=segmentKernel_->stressMoments(helper,piece,center,width);
                        for(int j=0;j<P;++j)
                        {
                            for(int c=0;c<9;++c)
                            {
                                h(b*P+j,c)+=values(j,c);
                            }
                        }
                    }
                }
                Kokkos::deep_copy(ex,coefficients,h);interaction->setIntegratedLeafMoments(coefficients);
                std::cout<<"  FMM stress: segment-integrated coefficients, "<<effective.size()<<" occupied leaves, polynomial Gauss order="<<helper.stressRuleSize()<<std::endl;
            }

        }
        DislocationStressArrays<3> evaluate()
        {
            if(options_.segmentShortRange) {
      interaction->evaluateLong(strengths,output);
    } else {
      interaction->evaluate(strengths,output);
    }
    auto host=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),output);
            leafLocals_=interaction->hostLeafLocals();
            DislocationStressArrays<3> result; result.reset(originalTargets_);
            if(options_.segmentShortRange){
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for(std::size_t t=0;t<originalTargets_;++t){Eigen::Vector3d x;for(int c=0;c<3;++c) {
          x[c]=targetSnapshot_.positions.values[c][t];
        }
        const auto value=stressAt(x);for(int i=0;i<3;++i) {
          for(int j=0;j<3;++j) {
            result.stresses.values[3*i+j][t]=value(i,j);
          }
        }
      }
                return result;
            }
            for(std::size_t t=0;t<originalTargets_;++t){int c=0;for(int i=0;i<3;++i) {
        for(int j=i;j<3;++j){result.stresses.values[3*i+j][t]=result.stresses.values[3*j+i][t]=host(t,c++);}
      }
    }
            return result;
        }
        Eigen::Matrix3d stressAt(const Eigen::Vector3d& x) const
        {
            if(!x.allFinite()||tree_->revision()!=revision_||leafLocals_.empty()) {
      throw std::logic_error("Missing or stale FMM stress cache");
    }
    double sum[6]={};const auto longKernel=family_.longKernel();const auto shortKernel=family_.shortKernel();
            const int n=1<<tree_->depth(),p=options_.interpolationOrder,P=p*p*p;
            const Eigen::Vector3d lo=Eigen::Vector3d::Zero(),width=Eigen::Vector3d::Constant(1.0/n);
            const Eigen::Vector3d fractional=tree_->fractional(x);
            // FEM face quadrature can round just outside an exact mesh bound.
            const double boundTolerance=boxBoundaryTolerance(tree_->lower(),tree_->upper());
            const bool inside=tree_->contains(x,boundTolerance);
            if(!inside && options_.screened && options_.boxShortRange) {
      throw std::out_of_range("Ewald stress target lies outside the prepared box tree; rebuild with bounds containing the target");
    }
    int cell[3]={};
            if(inside)
            {
                double w[3][6]={};
                for(int c=0;c<3;++c)
                {
                    cell[c]=std::clamp(int(fractional[c]/width[c]),0,n-1);
                    const double q=2*(fractional[c]-(cell[c]+.5)*width[c])/width[c];
                    for(int j=0;j<p;++j)
                    {
                        const double node=std::cos(3.14159265358979323846*j/(p-1));w[c][j]=1;
                        for(int k=0;k<p;++k) {
            if(j!=k){const double other=std::cos(3.14159265358979323846*k/(p-1));w[c][j]*=(q-other)/(node-other);}
          }
        }
                }
                const std::size_t leaf=cell[0]+n*(cell[1]+n*cell[2]);
                for(int j=0;j<P;++j){const double weight=w[0][j%p]*w[1][(j/p)%p]*w[2][j/(p*p)];for(int c=0;c<6;++c) {
          sum[c]+=weight*leafLocals_[(leaf*P+j)*6+c];
        }
      }
            }
            const auto& points=index_.sources();
            const bool stencil=inside && (!options_.screened || options_.boxShortRange);
            const std::size_t targetLeaf=cell[0]+n*(cell[1]+n*cell[2]);
            const auto& nearEntries=nearNeighborhood_->rows[targetLeaf];
            const std::size_t leafCount=std::size_t(n)*n*n;
            const std::size_t candidateCount=stencil?nearEntries.size():leafCount*options_.localImageShifts.size();
            for(std::size_t candidate=0;candidate<candidateCount;++candidate)
            {
                const auto leaf=stencil?nearEntries[candidate].source:candidate%leafCount;
                const auto image=stencil?nearEntries[candidate].image:candidate/leafCount;
                const Eigen::Vector3d shiftedTarget=x-options_.localImageShifts[image];
                const bool adjacent=inside&&!tree_->admissible(tree_->cellCenter(tree_->depth(),targetLeaf)-tree_->cellCenter(tree_->depth(),leaf)-options_.localImageShifts[image],tree_->depth());
                const bool near=!inside||(!(options_.screened && options_.interpolateScreenedNear)&&adjacent);
                const double distance=std::max(0.0,(shiftedTarget-tree_->cellCenter(tree_->depth(),leaf)).norm()-.5*tree_->cellDiameter(tree_->depth()));
                const double dist2=distance*distance;
                const bool local=options_.screened && (options_.boxShortRange?adjacent:dist2<=options_.shortCutoff*options_.shortCutoff);
                if(!near&&!local) {
            continue;
          }
          if(options_.segmentShortRange && (local||near)){Eigen::Matrix3d value=Eigen::Matrix3d::Zero();for(const auto& piece:segmentBins_[leaf]) {
              value+=segmentKernel_->stressNear(shiftedTarget,piece);
            }
            int c=0;for(int i=0;i<3;++i) {
              for(int j=i;j<3;++j) {
                sum[c++]+=value(i,j);
              }
            }
          }
                if(options_.segmentShortRange) {
            continue;
          }
          for(std::size_t k=points.offsets[lexToBox_[leaf]];k<points.offsets[lexToBox_[leaf]+1];++k)
                {
                    const auto s=points.permutation[k];double r[3],q[9];
                    for(int c=0;c<3;++c) {
              r[c]=shiftedTarget[c]-sources_.positions.values[c][s];
            }
            for(int c=0;c<9;++c) {
              q[c]=sources_.weightedBurgersTensors.values[c][s];
            }
            if(near){if(options_.screened) {
                longKernel.add(r,q,sum);
              } else {
                family_.fullKernel().add(r,q,sum);
              }
            }
                    if(!options_.segmentShortRange&&local&&(options_.boxShortRange||shortKernel.distance(r)<=options_.shortCutoff)) {
              shortKernel.add(r,q,sum);
            }
          }
            }
    Eigen::Matrix3d result;int c=0;for(int i=0;i<3;++i) {
      for(int j=i;j<3;++j) {
        result(i,j)=result(j,i)=sum[c++];
      }
    }
    return result;
        }
    };
    static ScreenedDislocationOptions concentrationSplitOptions(const BoxTree<3>& tree,
        const std::vector<ClimbKernelChannel>& channels,const ScreenedDislocationOptions& supplied)
    {
        auto options=supplied;
        if(!options.segmentSplits)
        {
            std::vector<SegmentSplitSnapshot::Endpoints> geometry;
            for(const auto& channel:channels)
            {
                for(const auto& segment:channel.segments)
                {
                    geometry.emplace_back(segment.p0,segment.p1);
                }
            }
            options.segmentSplits=std::make_shared<SegmentSplitSnapshot>(tree,SegmentSplitSnapshot::geometry(geometry));
        }
        if(!options.segmentSplits->matches(tree))
        {
            throw std::logic_error("Stale shared segment split tree");
        }
        return options;
    }

    struct DislocationConcentrationEvaluator::Impl
    {
        Registrations registrations;
        Ex ex;
        std::shared_ptr<const BoxTree<3>> tree;
        std::uint64_t revision;
        std::size_t nodes,species,targets;
        bool longOnly;

        struct SparseMap
        {
            Indices offsets,columns;
            View values;
            void upload(const Ex& ex,const Eigen::SparseMatrix<double,Eigen::RowMajor>& matrix)
            {
                offsets=Indices("segment map offsets",matrix.rows()+1,1);
                columns=Indices("segment map columns",matrix.nonZeros(),1);
                values=View("segment map values",matrix.nonZeros(),1);
                auto o=Kokkos::create_mirror(offsets),c=Kokkos::create_mirror(columns);
                auto w=Kokkos::create_mirror(values);
                for(Eigen::Index i=0;i<=matrix.rows();++i)
                {
                    o(i,0)=matrix.outerIndexPtr()[i];
                }
                for(Eigen::Index i=0;i<matrix.nonZeros();++i)
                {
                    c(i,0)=matrix.innerIndexPtr()[i];w(i,0)=matrix.valuePtr()[i];
                }
                Kokkos::deep_copy(ex,offsets,o);Kokkos::deep_copy(ex,columns,c);Kokkos::deep_copy(ex,values,w);
            }
            void apply(const Ex& ex,View velocity,View out,bool add=false) const
            {
                const auto o=offsets,c=columns;const auto w=values;
                Kokkos::parallel_for("segment concentration map",Kokkos::RangePolicy<Ex,Kokkos::IndexType<std::size_t>>(ex,0,out.extent(0)),KOKKOS_LAMBDA(std::size_t i)
                {
                    double value=0;
                    for(std::size_t j=o(i,0);j<o(i+1,0);++j)
                    {
                        value+=w(j,0)*velocity(c(j,0),0);
                    }
                    out(i,0)=value+(add?out(i,0):0);
                });
            }
        };
        struct Channel
        {
            std::size_t species;
            std::vector<std::size_t> rows;
            std::shared_ptr<KokkosFMMInteraction<HelmholtzDiffusionFamily>> interaction;
            Indices ids;
            View weights,strengths,output,moments,velocity;
            std::vector<std::size_t> nodeIDs;
            SparseMap momentMap,localMap;
            bool segmentBased=false,addLocal=false;
        };
        std::vector<Channel> channels;
        Impl(std::shared_ptr<FMMManager> manager,std::shared_ptr<const BoxTree<3>> t,
             std::size_t n,std::size_t ns,const std::vector<ClimbKernelChannel>& sources,
             const std::vector<ConcentrationTarget>& points,const ScreenedDislocationOptions& supplied,bool filter,bool onlyLong)
            :registrations(std::move(manager)),tree(std::move(t)),revision(tree?tree->revision():0),nodes(n),species(ns),targets(points.size()),longOnly(onlyLong)
        {
            if(!tree)
            {
                throw std::invalid_argument("Missing concentration tree");
            }
            const auto options=concentrationSplitOptions(*tree,sources,supplied);
            if(options.screened && options.coreType==CoreType::Cai)
            {
                throw std::invalid_argument("ewaldFMM requires coreType=Helmholtz");
            }
            if(!tree||species==0) {
      throw std::invalid_argument("Invalid concentration evaluator geometry");
    }
    for(std::size_t ci=0;ci<sources.size();++ci)
            {
                auto source=sources[ci];
                if(source.species>=species) {
        throw std::invalid_argument("Invalid concentration species");
      }
      Channel ch;ch.species=source.species;Positions targetPositions;
                for(std::size_t i=0;i<points.size();++i)
                {
                    if(filter&&std::find(points[i].grains.begin(),points[i].grains.end(),source.grain)==points[i].grains.end()) {
          continue;
        }
        ch.rows.push_back(i);for(int c=0;c<3;++c) {
          targetPositions[c].push_back(points[i].position[c]);
        }
      }
                if(!source.sources.empty()&&!source.segments.empty())
                {
                    throw std::invalid_argument("Mixed particle and segment concentration sources");
                }
                if(ch.rows.empty()||(source.sources.empty()&&source.segments.empty())) {
        continue;
      }
                std::set<std::size_t> nodeIDs;
                for(const auto& entry:source.sources) { nodeIDs.insert(entry.nodes.begin(),entry.nodes.end()); }
                for(const auto& entry:source.segments) { nodeIDs.insert(entry.nodes.begin(),entry.nodes.end()); }
                ch.nodeIDs.assign(nodeIDs.begin(),nodeIDs.end());
                std::map<std::size_t,std::size_t> localIDs;
                for(std::size_t i=0;i<ch.nodeIDs.size();++i)
                {
                    if(ch.nodeIDs[i]>=nodes) { throw std::invalid_argument("Invalid concentration source node"); }
                    localIDs.emplace(ch.nodeIDs[i],i);
                }
                const auto remap=[&](auto& entry)
                {
                    for(auto& node:entry.nodes) { node=localIDs.at(node); }
                };
                for(auto& entry:source.sources) { remap(entry); }
                for(auto& entry:source.segments) { remap(entry); }
                ch.velocity=View("grain concentration velocities",ch.nodeIDs.size(),1);
                ch.segmentBased=!source.segments.empty();
                ch.addLocal=!onlyLong||!options.screened;
                std::unique_ptr<SegmentConcentrationPlan> segmentPlan;
                Positions src;
                auto opt=interactionOptions(options);
                opt.buildShortNeighbors=!onlyLong;
                if(ch.segmentBased)
                {
                    segmentPlan=std::make_unique<SegmentConcentrationPlan>(*tree,source,ch.nodeIDs.size(),targetPositions,options,ch.addLocal);
                    src=segmentPlan->markers;
                    if(profilingEnabled())
        {
            std::cout<<"  FMM concentration: "<<source.segments.size()<<" source segments, "
                             <<src[0].size()<<" occupied leaves, segment-integrated source coefficients"<<std::endl;
        }
                    opt.buildShortNeighbors=false;
                    opt.externalNear=!options.screened;
                    ch.momentMap.upload(ex,segmentPlan->moments);
                    ch.moments=View("segment concentration moments",segmentPlan->moments.rows(),1);
                    if(ch.addLocal)
                    {
                        ch.localMap.upload(ex,segmentPlan->local);
                    }
                }
                else
                {
                    src=positions(source.sources,ch.nodeIDs.size());
                    uploadMap(ex,source.sources,ch.ids,ch.weights);
                }
                HelmholtzDiffusionFamily family(source.diffusion,options.coreWidth,options.screenWidth,options.coreType,options.screened);
                auto handle=registrations.addKernel("concentration/"+std::to_string(ci),family);
                ch.interaction=registrations.addInteraction<KokkosFMMInteraction<HelmholtzDiffusionFamily>>("concentration/interaction/"+std::to_string(ci),handle,ex,tree,src,targetPositions,opt);
                ch.strengths=View("concentration strengths",src[0].size(),1);
                ch.output=View("concentration outputs",ch.rows.size(),1);
                channels.push_back(std::move(ch));
            }
        }
        Eigen::VectorXd evaluate(std::size_t speciesID,const Eigen::VectorXd& v)
        {
            // Opt-in because fences perturb asynchronous device execution.
            const bool profile=longOnly && profilingEnabled();
            using Clock=std::chrono::steady_clock;
            auto stamp=Clock::now();
            const auto start=stamp;
            std::array<double,5> times{}; // upload, source map, long plan, local map, download
            std::array<double,4> planTimes{};
            std::array<double,3> translationTimes{};
            const auto checkpoint=[&](std::size_t bucket)
            {
                if(profile)
                {
                    ex.fence();
                    const auto now=Clock::now();
                    times[bucket]+=std::chrono::duration<double>(now-stamp).count();
                    stamp=now;
                }
            };
            if(speciesID>=species||v.size()!=Eigen::Index(nodes)||!v.allFinite()) {
      throw std::invalid_argument("Invalid concentration velocities");
    }
    if(tree->revision()!=revision) {
      throw std::logic_error("Stale concentration geometry");
    }
            Eigen::VectorXd result=Eigen::VectorXd::Zero(targets);
            for(auto& ch:channels) {
      if(ch.species==speciesID)
            {
                auto host=Kokkos::create_mirror(ch.velocity);
                for(std::size_t i=0;i<ch.nodeIDs.size();++i) { host(i,0)=v(ch.nodeIDs[i]); }
                Kokkos::deep_copy(ex,ch.velocity,host);
                checkpoint(0);
                ch.interaction->setProfiling(profile);
                if(ch.segmentBased)
                {
                    ch.momentMap.apply(ex,ch.velocity,ch.moments);
                    ch.interaction->setIntegratedLeafMoments(ch.moments);
                    checkpoint(1);
                    ch.interaction->evaluateLong(ch.strengths,ch.output);
                    checkpoint(2);
                    if(ch.addLocal)
                    {
                        ch.localMap.apply(ex,ch.velocity,ch.output,true);
                    }
                }
                else
                {
                const auto ids=ch.ids;const auto weights=ch.weights,q=ch.strengths,vv=ch.velocity;
                Kokkos::parallel_for("concentration source map",Kokkos::RangePolicy<Ex,Kokkos::IndexType<std::size_t>>(ex,0,q.extent(0)),KOKKOS_LAMBDA(const std::size_t i)
                {q(i,0)=weights(i,0)*vv(ids(i,0),0)+weights(i,1)*vv(ids(i,1),0);});
                checkpoint(1);
                if(longOnly) {
          ch.interaction->evaluateLong(q,ch.output);
        } else {
          ch.interaction->evaluate(q,ch.output);
        }
                checkpoint(2);
                }
                checkpoint(3);
                if(profile)
                {
                    for(std::size_t j=0;j<planTimes.size();++j)
                    {
                        planTimes[j]+=0.001*ch.interaction->lastTimings()[j];
                    }
                    for(std::size_t j=0;j<translationTimes.size();++j)
                    {
                        translationTimes[j]+=0.001*ch.interaction->translationTimings()[j];
                    }
                }
        auto values=Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),ch.output);
                for(std::size_t i=0;i<ch.rows.size();++i) {
          result(ch.rows[i])+=values(i,0);
        }
                checkpoint(4);
      }
    }
    if(profile)
    {
        std::ostringstream message;
        message<<std::scientific<<std::setprecision(9)
               <<"\n  FMM operator detail: upload="<<times[0]<<", source_map="<<times[1]
               <<", plan="<<times[2]<<", local_map="<<times[3]<<", download="<<times[4]
               <<", cache="<<planTimes[0]<<", upward="<<planTimes[1]
               <<", m2l="<<translationTimes[0]<<", l2l="<<translationTimes[1]
               <<", near_m2l="<<translationTimes[2]<<", target="<<planTimes[3]
               <<", total="<<std::chrono::duration<double>(Clock::now()-start).count();
        std::cout<<message.str();
    }
    return result;
        }
    };
    struct ScreenedClimbAssembly::Impl
    {
        Registrations registrations;
        Ex ex;
        std::shared_ptr<const BoxTree<3>> tree;
        std::uint64_t revision;
        std::size_t nodes,species;
        bool screened;
        std::vector<Eigen::SparseMatrix<double>> local;
        std::vector<Eigen::VectorXd> rhs;
        std::unique_ptr<DislocationConcentrationEvaluator> evaluator;
        std::vector<ClimbQuadratureMap> targetMaps;
        std::vector<std::vector<std::pair<std::size_t,double>>> targetRows;
        Impl(std::shared_ptr<FMMManager> m,std::shared_ptr<const BoxTree<3>> t,
             const ClimbAssemblySnapshot& snapshot,const ScreenedDislocationOptions& supplied)
            :registrations(std::move(m)),tree(std::move(t)),revision(tree?tree->revision():0),nodes(snapshot.nodeCount),species(snapshot.speciesCount),screened(supplied.screened),rhs(snapshot.rhs),targetMaps(snapshot.targets)
        {
            if(!tree)
            {
                throw std::invalid_argument("Missing climb tree");
            }
            const auto opt=concentrationSplitOptions(*tree,snapshot.channels,supplied);
            if(!std::isfinite(opt.coreWidth)||!std::isfinite(opt.screenWidth)||opt.coreWidth<=0||
               opt.screenWidth<=0||opt.coreWidth/opt.screenWidth>.25||
               ((!opt.screened||!opt.boxShortRange)&&(!std::isfinite(opt.shortCutoff)||opt.shortCutoff<=0))||opt.interpolationOrder<2||opt.interpolationOrder>6) {
      throw std::invalid_argument("Invalid screened climb settings");
    }
    if(!tree || species==0 || rhs.size()!=species || nodes>std::size_t(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("Invalid climb snapshot dimensions");
    }
    for(const auto& f:rhs) {
      if(f.size()!=Eigen::Index(nodes)||!f.allFinite()) {
        throw std::invalid_argument("Invalid climb RHS");
      }
    }
    const auto targets=positions(snapshot.targets,nodes);
            targetRows.resize(nodes);
            for(std::size_t i=0;i<targetMaps.size();++i)
            {
                for(int k=0;k<2;++k)
                {
                    targetRows.at(targetMaps[i].nodes[k]).emplace_back(i,targetMaps[i].weights[k]);
                }
            }
            std::vector<ConcentrationTarget> fieldPoints;
            for(const auto& point:snapshot.targets){ConcentrationTarget target;target.position=point.position;target.grains={point.grain};fieldPoints.push_back(target);}
            SolverStageTimer setupStage("Climb setup: FMM plans",true);
            evaluator=std::make_unique<DislocationConcentrationEvaluator>(registrations.manager,tree,nodes,species,snapshot.channels,fieldPoints,opt,true,true);
            setupStage.start("Climb setup: local matrix construction");
            std::vector<std::vector<Eigen::Triplet<double>>> entries(species);
            for(std::size_t ci=0;ci<snapshot.channels.size();++ci)
            {
                const auto& input=snapshot.channels[ci];
                if(input.species>=species) {
        throw std::invalid_argument("Invalid climb species");
      }
                if(!input.segments.empty())
                {
                    Positions grainTargets;
                    std::vector<std::size_t> grainRows;
                    for(std::size_t i=0;i<snapshot.targets.size();++i)
                    {
                        if(snapshot.targets[i].grain==input.grain)
                        {
                            grainRows.push_back(i);
                            for(int c=0;c<3;++c) { grainTargets[c].push_back(snapshot.targets[i].position[c]); }
                        }
                    }
                    SegmentConcentrationPlan plan(*tree,input,nodes,grainTargets,opt,true,false);
                    auto& triplets=entries[input.species];
                    std::vector<std::vector<Eigen::Triplet<double>>> rows(grainRows.size());
                    independentParallelFor(rows.size(),[&](std::size_t ti)
                    {
                        const auto& target=snapshot.targets[grainRows[ti]];
                        for(Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(plan.local,ti);it;++it)
                        {
                            for(int k=0;k<2;++k)
                            {
                                rows[ti].emplace_back(target.nodes[k],it.col(),target.weights[k]*it.value());
                            }
                        }
                    });
                    appendOrderedRows(triplets,rows);
                    continue;
                }
      const auto sources=positions(input.sources,nodes);
                HelmholtzDiffusionFamily family(input.diffusion,opt.coreWidth,opt.screenWidth,opt.coreType,opt.screened);
                // Assemble only geometrically local point pairs, scatter four endpoint terms.
                std::unique_ptr<ShortRangeParticleNeighbors> neighbors;
                std::unique_ptr<LeafBoxNeighbors> boxNeighbors;
                if((opt.screened && opt.boxShortRange)||opt.localImageShifts.size()!=1||!opt.localImageShifts.front().isZero()) {
        boxNeighbors=std::make_unique<LeafBoxNeighbors>(*tree,sources,targets,opt.localImageShifts);
      } else {
        neighbors=std::make_unique<ShortRangeParticleNeighbors>(sources,targets,family.searchRadius(opt.shortCutoff));
      }
      const auto kernel=family.shortKernel();auto& triplets=entries[input.species];
                for(std::size_t ti=0;ti<snapshot.targets.size();++ti)
                {
                    const auto& target=snapshot.targets[ti];
                    if(target.grain!=input.grain) { continue; }
                    const auto accumulate=[&](std::size_t sourceID,std::size_t image)
                    {
                        const auto& source=input.sources[sourceID];
                        double r[3];for(int c=0;c<3;++c) {
            r[c]=target.position[c]-source.position[c]-opt.localImageShifts[image][c];
          }
          if(!boxNeighbors && kernel.distance(r)>opt.shortCutoff) {
            return;
          }
          const double g=opt.screened?kernel.value(r):family.fullKernel().value(r);
                        for(int i=0;i<2;++i) {
            for(int k=0;k<2;++k)
                        {
                            const double value=target.weights[i]*g*source.weights[k];
                            if(value!=0) {
                triplets.emplace_back(target.nodes[i],source.nodes[k],value);
              }
            }
          }
        };
                    if(boxNeighbors) {
          boxNeighbors->forEachImage(ti,accumulate);
        } else {
          for(std::size_t j=neighbors->offsets[ti];j<neighbors->offsets[ti+1];++j) {
            accumulate(neighbors->pairs[j].source,0);
          }
        }
      }
            }
            for(auto& triplets:entries){Eigen::SparseMatrix<double> matrix(nodes,nodes);matrix.setFromTriplets(triplets.begin(),triplets.end());matrix.makeCompressed();local.push_back(std::move(matrix));}
        }
        void check(std::size_t s,const Eigen::VectorXd& v) const
        {
            if(s>=species||v.size()!=Eigen::Index(nodes)||!v.allFinite()) {
      throw std::invalid_argument("Invalid climb operator input");
    }
    if(tree->revision()!=revision) {
      throw std::logic_error("Climb tree changed; reassemble geometry");
    }
  }
        Eigen::VectorXd applyLong(std::size_t s,const Eigen::VectorXd& v)
        {
            check(s,v);
            const auto concentration=evaluator->evaluate(s,v);
            Eigen::VectorXd out=Eigen::VectorXd::Zero(nodes);
            independentParallelFor(nodes,[&](std::size_t node)
            {
                double value=0.0;
                for(const auto& term:targetRows[node]) { value+=term.second*concentration(term.first); }
                out(node)=value;
            });
    return out;
        }
    };
#else
    struct FMMStressEvaluator::Impl {};
    struct ScreenedClimbAssembly::Impl {};
    struct DislocationConcentrationEvaluator::Impl {};
#endif
    DislocationConcentrationEvaluator::DislocationConcentrationEvaluator(std::shared_ptr<FMMManager> m,std::shared_ptr<const BoxTree<3>> t,
        std::size_t nodes,std::size_t species,const std::vector<ClimbKernelChannel>& sources,
        const std::vector<ConcentrationTarget>& targets,const ScreenedDislocationOptions& options,bool filter,bool longOnly)
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        impl_=std::make_unique<Impl>(std::move(m),std::move(t),nodes,species,sources,targets,options,filter,longOnly);
#else
        throw std::logic_error("Concentration FMM requires Kokkos");
#endif
    }
    DislocationConcentrationEvaluator::~DislocationConcentrationEvaluator()=default;
    Eigen::VectorXd DislocationConcentrationEvaluator::evaluate(std::size_t species,const Eigen::VectorXd& v)
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->evaluate(species,v);
#else
        throw std::logic_error("Concentration FMM requires Kokkos");
#endif
    }
    FMMStressEvaluator::FMMStressEvaluator(std::shared_ptr<FMMManager> m,std::shared_ptr<const BoxTree<3>> t,
        const DislocationSourceArrays<3>& s,const DislocationTargetArrays<3>& f,double mu,double nu,const ScreenedDislocationOptions& o,const std::vector<helmholtzSegment::Segment>& segments)
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        if((o.coreType==CoreType::Cai))
        {
            if(o.screened) {
      throw std::invalid_argument("Cai stress requires multipoleMode=FMM (no Cai Ewald split implemented)");
    }
    impl_=std::make_unique<StressImpl<CaiStressKernelFamily>>(std::move(m),std::move(t),s,f,mu,nu,o,segments);
        }
        else {
    impl_=std::make_unique<StressImpl<HelmholtzStressKernelFamily>>(std::move(m),std::move(t),s,f,mu,nu,o,segments);
  }
#else
        throw std::logic_error("Screened stress requires Kokkos");
#endif
    }
    FMMStressEvaluator::~FMMStressEvaluator()=default;
    DislocationStressArrays<3> FMMStressEvaluator::evaluate()
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->evaluate();
#else
        throw std::logic_error("Screened stress requires Kokkos");
#endif
    }
    Eigen::Matrix3d FMMStressEvaluator::stressAt(const Eigen::Vector3d& x) const
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->stressAt(x);
#else
        throw std::logic_error("Screened stress requires Kokkos");
#endif
    }
    ScreenedClimbAssembly::ScreenedClimbAssembly(std::shared_ptr<FMMManager> m,std::shared_ptr<const BoxTree<3>> t,
        const ClimbAssemblySnapshot& s,const ScreenedDislocationOptions& o)
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        impl_=std::make_unique<Impl>(std::move(m),std::move(t),s,o);
#else
        throw std::logic_error("Screened climb requires Kokkos");
#endif
    }
    ScreenedClimbAssembly::~ScreenedClimbAssembly()=default;
    std::size_t ScreenedClimbAssembly::nodeCount() const
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->nodes;
#else
        return 0;
#endif
    }
    std::size_t ScreenedClimbAssembly::speciesCount() const
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->species;
#else
        return 0;
#endif
    }
    const Eigen::SparseMatrix<double>& ScreenedClimbAssembly::shortMatrix(std::size_t s) const
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        if(impl_->tree->revision()!=impl_->revision) {
    throw std::logic_error("Climb tree changed; reassemble geometry");
  }
  return impl_->local.at(s);
#else
        throw std::logic_error("Screened climb requires Kokkos");
#endif
    }
    const Eigen::VectorXd& ScreenedClimbAssembly::rhs(std::size_t s) const
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->rhs.at(s);
#else
        throw std::logic_error("Screened climb requires Kokkos");
#endif
    }
    Eigen::VectorXd ScreenedClimbAssembly::applyLong(std::size_t s,const Eigen::VectorXd& v)
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return impl_->applyLong(s,v);
#else
        throw std::logic_error("Screened climb requires Kokkos");
#endif
    }
    Eigen::VectorXd ScreenedClimbAssembly::apply(std::size_t s,const Eigen::VectorXd& v)
    {auto out=applyLong(s,v);
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        if(impl_->screened) {
    out+=shortMatrix(s)*v;
  }
#endif
        return out;}
}
