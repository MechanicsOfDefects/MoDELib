/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_DislocationFieldEvaluator_H_
#define model_DislocationFieldEvaluator_H_
#include <DislocationDynamics/DislocationFieldBuffers.h>
#include <FMM/BoxParticleIndex.h>
#include <DislocationDynamics/ScreenedDislocationAssembly.h>
#include <memory>
namespace model
{
    template<int dim> class DislocationNetwork;
    // Network-owned, explicit snapshots. Rebuild after any geometry/topology
    // change before evaluating fields. Screened stress uses an explicit prepared plan.
    template<int dim>
    class DislocationFieldEvaluator
    {
        struct State;
        struct Device;
        const DislocationNetwork<dim>& network_;
        const std::shared_ptr<BoxTree<dim>> tree_;
        std::unique_ptr<State> state_;
        std::unique_ptr<Device> device_;
    public:
        struct Timings { double samplingMilliseconds, assignmentMilliseconds; };
        DislocationFieldEvaluator(const DislocationNetwork<dim>&,std::shared_ptr<BoxTree<dim>>);
        ~DislocationFieldEvaluator();
        DislocationFieldEvaluator(const DislocationFieldEvaluator&)=delete;
        DislocationFieldEvaluator& operator=(const DislocationFieldEvaluator&)=delete;
        Timings rebuild(double sourceDensity,double targetDensity,DislocationTargetPurpose);
        bool hasSnapshot() const;
        const DislocationFieldBuffers<dim>& buffers() const;
        const BoxParticleIndex<dim>& leafIndex() const;
        const BoxTree<dim>& tree() const;
        static bool kokkosEnabled();
        // Lazy device allocation. Caller owns Kokkos lifetime: releaseDevice()
        // (or destroy the network) before finalizing Kokkos.
        void upload();
        // Free-space Lazar stress at the snapshot's independently sampled target points.
        void prepareFMMStress(double mu,double nu,const ScreenedDislocationOptions&);
        DislocationStressArrays<dim> evaluateFMMStress();
        Eigen::Matrix<double,dim,dim> fmmStressAt(const Eigen::Matrix<double,dim,1>&) const;
        bool hasFMMStress() const;
        void releaseDevice();
        void clear();
    };
    extern template class DislocationFieldEvaluator<3>;
}
#endif
