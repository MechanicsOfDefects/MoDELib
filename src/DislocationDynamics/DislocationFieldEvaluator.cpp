#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
#include <DislocationDynamics/KokkosDislocationFieldBuffers.h>
#include <FMM/KokkosBoxParticleIndex.h>
#endif
#include <DislocationDynamics/DislocationFieldEvaluator.h>
#include <DislocationDynamics/DislocationFieldSampling.h>
#include <DislocationDynamics/DislocationNetwork.h>
#include <chrono>
namespace model
{
    template<int dim> struct DislocationFieldEvaluator<dim>::State
    {
        DislocationFieldBuffers<dim> buffers;
        BoxParticleIndex<dim> index;
        std::unique_ptr<FMMStressEvaluator> stress;
    };
    template<int dim> struct DislocationFieldEvaluator<dim>::Device
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        KokkosDislocationFieldBuffers<dim> buffers;
        KokkosBoxParticleIndex<dim> index;
#endif
    };
    template<int dim>
    DislocationFieldEvaluator<dim>::DislocationFieldEvaluator(const DislocationNetwork<dim>& network,
                                                             std::shared_ptr<BoxTree<dim>> tree)
        :network_(network),tree_(std::move(tree))
    {
        if(!tree_) { throw std::invalid_argument("Field evaluator requires a box tree"); }
    }
    template<int dim> DislocationFieldEvaluator<dim>::~DislocationFieldEvaluator()=default;
    template<int dim> typename DislocationFieldEvaluator<dim>::Timings
    DislocationFieldEvaluator<dim>::rebuild(double sourceDensity,double targetDensity,DislocationTargetPurpose purpose)
    {
        if(!std::isfinite(sourceDensity) || sourceDensity<=0 || !std::isfinite(targetDensity) || targetDensity<=0)
        { throw std::invalid_argument("Positive finite quadrature densities required"); }
        using Clock=std::chrono::steady_clock;
        auto next=std::make_unique<State>();
        const auto start=Clock::now();
        for(const auto& entry:network_.networkLinks())
        {
            const auto segment=entry.second.lock();
            if(!segment) { throw std::logic_error("Expired network segment during field sampling"); }
            DislocationFieldSampling<dim>::append(next->buffers,*segment,sourceDensity,targetDensity,purpose);
        }
        const auto sampled=Clock::now();
        next->index.rebuild(*tree_,next->buffers.sources().positions.values,next->buffers.targets().positions.values);
        const auto indexed=Clock::now();
        // Publish only after sampling and assignment both succeed.
        releaseDevice();
        state_=std::move(next);
        return {std::chrono::duration<double,std::milli>(sampled-start).count(),
                std::chrono::duration<double,std::milli>(indexed-sampled).count()};
    }
    template<int dim> bool DislocationFieldEvaluator<dim>::hasSnapshot() const { return bool(state_); }
    template<int dim> const DislocationFieldBuffers<dim>& DislocationFieldEvaluator<dim>::buffers() const
    {
        if(!state_) { throw std::logic_error("Rebuild field evaluator before accessing its snapshot"); }
        return state_->buffers;
    }
    template<int dim> const BoxParticleIndex<dim>& DislocationFieldEvaluator<dim>::leafIndex() const
    {
        if(!state_ || !state_->index.matches(*tree_)) { throw std::logic_error("Missing or stale field leaf index"); }
        return state_->index;
    }
    template<int dim> const BoxTree<dim>& DislocationFieldEvaluator<dim>::tree() const { return *tree_; }
    template<int dim> bool DislocationFieldEvaluator<dim>::kokkosEnabled()
    {
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        return true;
#else
        return false;
#endif
    }
    template<int dim> void DislocationFieldEvaluator<dim>::upload()
    {
        leafIndex();
#ifdef MODELIB_FIELD_EVALUATOR_KOKKOS
        if(!device_) { device_=std::make_unique<Device>(); }
        device_->buffers.uploadSources(buffers().sources());
        device_->buffers.uploadTargets(buffers().targets());
        device_->index.upload(leafIndex(),*tree_);
#else
        throw std::logic_error("MoDELib was built without Kokkos");
#endif
    }
    template<int dim>
    void DislocationFieldEvaluator<dim>::prepareFMMStress(double mu,double nu,const ScreenedDislocationOptions& options)
    {
        leafIndex();
        if(network_.ddBase.periodicLatticeBasis.cols()!=0)
        { throw std::logic_error("Screened stress periodic correction is not implemented"); }
        std::vector<helmholtzSegment::Segment> segments;
        if(options.segmentShortRange) {
    for(const auto& entry:network_.networkLinks())
        {const auto seg=entry.second.lock();if(!seg) {
        throw std::logic_error("Expired stress segment");
      }
      if(!seg->hasZeroBurgers()) {
        segments.push_back({seg->source->get_P(),seg->sink->get_P(),seg->burgers()});
      }
    }
  }
  auto sharedOptions=options;
        sharedOptions.segmentSplits=network_.segmentSplits();
        auto next=std::make_unique<FMMStressEvaluator>(network_.ddBase.fmmContext->manager,tree_,
            buffers().sources(),buffers().targets(),mu,nu,sharedOptions,segments);
        next->evaluate(); // Publish complete, read-only leaf locals for threaded field queries.
        state_->stress=std::move(next);
    }
    template<int dim> DislocationStressArrays<dim> DislocationFieldEvaluator<dim>::evaluateFMMStress()
    {
        leafIndex();
        if(!state_->stress) {
    throw std::logic_error("Prepare screened stress after rebuilding the snapshot");
  }
  return state_->stress->evaluate();
    }
    template<int dim> bool DislocationFieldEvaluator<dim>::hasFMMStress() const { return state_&&bool(state_->stress); }
    template<int dim> Eigen::Matrix<double,dim,dim> DislocationFieldEvaluator<dim>::fmmStressAt(const Eigen::Matrix<double,dim,1>& x) const
    {leafIndex();if(!hasFMMStress()) {
    throw std::logic_error("Prepare FMM stress before field queries");
  }
  return state_->stress->stressAt(x);}
    template<int dim> void DislocationFieldEvaluator<dim>::releaseDevice() { if(state_) {
    state_->stress.reset();
  }
  device_.reset(); }
    template<int dim> void DislocationFieldEvaluator<dim>::clear() { releaseDevice(); state_.reset(); }
    template class DislocationFieldEvaluator<3>;
}
