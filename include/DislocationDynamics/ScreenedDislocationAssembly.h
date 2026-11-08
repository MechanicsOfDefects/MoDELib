#include <DislocationDynamics/CoreType.h>
/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_ScreenedDislocationAssembly_H_
#define model_ScreenedDislocationAssembly_H_
#include <FMM/FMMManager.h>
#include <FMM/BoxTree.h>
#include <DislocationDynamics/DislocationFieldBuffers.h>
#include <Eigen/Sparse>
#include <memory>
#include <DislocationDynamics/SegmentSplitCache.h>
#include <DislocationDynamics/HelmholtzSegmentGeometry.h>
namespace model
{
    // Shared free-space settings for stress, concentration and climb operators.
    struct ScreenedDislocationOptions
    {
        std::shared_ptr<const SegmentSplitSnapshot> segmentSplits;
        bool boxShortRange=true;
        bool segmentShortRange=false; // Set by DD for ewaldFMM; particle-only tests remain available.
        int segmentChebyshevOrder=8,segmentTableResolution=64;
        bool interpolateScreenedNear=true; // Smooth long kernel uses self/neighbor M2L.
        CoreType coreType=CoreType::Helmholtz;
        bool screened=true; // False: full core-regularized kernel, no additive short correction.
        double coreWidth=1,screenWidth=4;
        double shortCutoff=32; // Plain-FMM preconditioner / radial reference path only.
        int interpolationOrder=4;
        std::size_t translationCacheBytes=256u*1024u*1024u;
    };
    class FMMStressEvaluator
    {
        struct Impl;
        template<class Family> struct StressImpl;
        std::unique_ptr<Impl> impl_;
    public:
        FMMStressEvaluator(std::shared_ptr<FMMManager>,std::shared_ptr<const BoxTree<3>>,
            const DislocationSourceArrays<3>&,const DislocationTargetArrays<3>&,
            double mu,double nu,const ScreenedDislocationOptions&,
            const std::vector<helmholtzSegment::Segment>& segments={});
        ~FMMStressEvaluator();
        FMMStressEvaluator(const FMMStressEvaluator&)=delete;
        FMMStressEvaluator& operator=(const FMMStressEvaluator&)=delete;
        DislocationStressArrays<3> evaluate();
        Eigen::Matrix3d stressAt(const Eigen::Vector3d&) const;
    };
    struct ClimbQuadratureMap
    {
        std::array<double,3> position{};
        std::array<std::size_t,2> nodes{};
        std::array<double,2> weights{}; // Shape * b cross tangent dot climbDirection * dL; source also /bias.
    };
    struct ClimbSourceSegment
    {
        Eigen::Vector3d p0,p1;
        std::array<std::size_t,2> nodes{};
        std::array<double,2> weights{}; // Endpoint b cross tangent dot climbDirection / bias, before shape and dL.
    };
    struct ClimbKernelChannel
    {
        std::size_t species=0,grain=0;
        Eigen::Matrix3d diffusion=Eigen::Matrix3d::Identity();
        std::vector<ClimbQuadratureMap> sources; // Particle reference path.
        std::vector<ClimbSourceSegment> segments; // Mutually exclusive with particle sources.
    };
    struct ConcentrationTarget
    {
        std::array<double,3> position{};
        std::vector<std::size_t> grains; // Used only when filterGrains=true.
    };
    // Reusable velocity -> source strengths -> concentration evaluator.
    class DislocationConcentrationEvaluator
    {
        struct Impl;
        std::unique_ptr<Impl> impl_;
    public:
        DislocationConcentrationEvaluator(std::shared_ptr<FMMManager>,std::shared_ptr<const BoxTree<3>>,
            std::size_t nodes,std::size_t species,const std::vector<ClimbKernelChannel>&,
            const std::vector<ConcentrationTarget>&,const ScreenedDislocationOptions&,
            bool filterGrains=false,bool longOnly=false);
        ~DislocationConcentrationEvaluator();
        Eigen::VectorXd evaluate(std::size_t species,const Eigen::VectorXd& velocities);
    };
    struct ClimbAssemblySnapshot
    {
        std::size_t nodeCount=0,speciesCount=0;
        std::vector<ClimbQuadratureMap> targets;
        std::vector<ClimbKernelChannel> channels; // Separate source-grain tensor per species.
        std::vector<Eigen::VectorXd> rhs;
    };
    // Assembles K_short and maps S,B; K_long is applied as B*FMM*S.
    // Does not solve, eliminate periodic DOFs, or retain live network pointers.
    class ScreenedClimbAssembly
    {
        struct Impl;
        std::unique_ptr<Impl> impl_;
    public:
        ScreenedClimbAssembly(std::shared_ptr<FMMManager>,std::shared_ptr<const BoxTree<3>>,
                             const ClimbAssemblySnapshot&,const ScreenedDislocationOptions&);
        ~ScreenedClimbAssembly();
        ScreenedClimbAssembly(const ScreenedClimbAssembly&)=delete;
        ScreenedClimbAssembly& operator=(const ScreenedClimbAssembly&)=delete;
        std::size_t nodeCount() const;
        std::size_t speciesCount() const;
        // In unscreened mode this is the local full-kernel preconditioner only.
        const Eigen::SparseMatrix<double>& shortMatrix(std::size_t species) const;
        const Eigen::VectorXd& rhs(std::size_t species) const;
        Eigen::VectorXd applyLong(std::size_t species,const Eigen::VectorXd& velocity);
        Eigen::VectorXd apply(std::size_t species,const Eigen::VectorXd& velocity);
    };
}
#endif
