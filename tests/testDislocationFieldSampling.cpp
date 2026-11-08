#ifdef MODELIB_TEST_EVALUATOR_KOKKOS
#include <cfloat>
#include <Kokkos_Core.hpp>
#endif
#include "ClusterPeriodicFixture.h"
#include <DislocationDynamics/DislocationFieldSampling.h>
#include <DislocationDynamics/GalerkinClimbSolver.h>
#ifdef MODELIB_TEST_EVALUATOR_KOKKOS
void checkEwaldInput()
{
    const auto folder=fs::temp_directory_path()/("ewald-input-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder);
    const auto dd=folder/"DD.txt",dc=folder/"DefectiveCrystal.txt";
    const auto writeFactor=[&](const std::string& factor)
    {std::ofstream(dd)<<"coreSize=2;\newaldWidthFactor="<<factor<<";\n";};
    std::ofstream(dc)<<"multipoleMode=ewaldFMM;\nfmmTreeDepth=2;\n";
    FMMContext<3> context(Eigen::Vector3d::Zero(),Eigen::Vector3d(800,1600,2400),dc.string());
    const auto boxSize=(context.boxTree->upper()-context.boxTree->lower()).minCoeff()/(1u<<context.fmmTreeDepth);
    writeFactor("0.125");
    require(DislocationMultipoleSettings(dd.string(),context.multipoleMode,boxSize).ewaldWidth==25,"Incorrect Ewald width from smallest leaf dimension");
    writeFactor("0.25");
    require(DislocationMultipoleSettings(dd.string(),context.multipoleMode,boxSize).ewaldWidth==50,"Factor did not scale width");
    require(context.boxTree->depth()==2,"Ewald factor changed user-selected tree depth");
    std::ofstream(dc)<<"multipoleMode=ewaldFMM;\nfmmTreeDepth=3;\n";
    FMMContext<3> finer(Eigen::Vector3d::Zero(),Eigen::Vector3d(800,1600,2400),dc.string());
    const double finerSize=(finer.boxTree->upper()-finer.boxTree->lower()).minCoeff()/(1u<<finer.fmmTreeDepth);
    require(DislocationMultipoleSettings(dd.string(),finer.multipoleMode,finerSize).ewaldWidth==25,"Depth did not scale Ewald width");
    for(const auto* value:{"0","-1","1","1.1","nan","garbage","0.01"})
    {
        writeFactor(value);bool rejected=false;
        try{DislocationMultipoleSettings bad(dd.string(),MultipoleMode::ewaldFMM,boxSize);}catch(const std::invalid_argument&){rejected=true;}
        require(rejected,"Invalid Ewald factor or derived core-width ratio accepted");
        // This parameter must not be parsed in either unscreened mode.
        require(DislocationMultipoleSettings(dd.string(),MultipoleMode::FMM,boxSize).ewaldWidth==8,"Plain FMM read Ewald-only factor");
        DislocationMultipoleSettings unused(dd.string(),MultipoleMode::none,boxSize);
    }
    std::ofstream(dd)<<"coreSize=2;\n";
    require(DislocationMultipoleSettings(dd.string(),MultipoleMode::ewaldFMM,boxSize).ewaldWidth==25,"Ewald default factor changed");
    std::ofstream(dd)<<"coreSize=2;\nfmmEwaldWidth=8;\n";
    bool retired=false;try{DislocationMultipoleSettings old(dd.string(),MultipoleMode::ewaldFMM,boxSize);}catch(const std::invalid_argument&){retired=true;}
    require(retired,"Retired absolute width was silently accepted");
    for(const auto& name:{std::string("Cai"),std::string("Helmholtz")})
    {
        std::ofstream(dd)<<"coreSize=2;\ncoreType="<<name<<";\newaldSegmentChebyshevOrder=10;\newaldSegmentTableResolution=32;\n";
        const auto settings=DislocationMultipoleSettings(dd.string(),MultipoleMode::none,boxSize);
        require((settings.coreType==CoreType::Cai)==(name=="Cai"),"none ignored coreType");
        require(settings.ewaldSegmentChebyshevOrder==10&&settings.ewaldSegmentTableResolution==32,"none ignored endpoint table settings");
    }
    for(const auto& key:{std::string("coreType=bad"),std::string("fmmStressCore=Cai")})
    {
        std::ofstream(dd)<<"coreSize=2;\n"<<key<<";\n";
        bool bad=false;
        try
        {
            DislocationMultipoleSettings settings(dd.string(),MultipoleMode::none,boxSize);
        }
        catch(const std::invalid_argument&)
        {
            bad=true;
        }
        require(bad,"Invalid or retired core selector accepted");
    }
    fs::remove_all(folder);
}
#endif
void check(const fs::path& root, int rank, int master=0)
{
    const auto folder = fixture(root, rank, 0, false);
    auto dc = read(folder / "inputFiles/DefectiveCrystal.txt");
    parameter(dc, "physics", "ClusterDynamics DislocationDynamics");
    // All-box long interpolation requires an accuracy check; use order six
    // for this strict direct-reference test, without changing the Ewald width.
    if(master==1) {
    dc+="\nfmmInterpolationOrder=6;\n";
  }
  dc += std::string("\nmultipoleMode=")+(master==2?"FMM":(master?"ewaldFMM":"none"))+";\n";
    if(rank>0)
    {
        parameter(dc,"periodicImageSize",rank==1?"1":"1 0 0");
    }
    std::ofstream(folder / "inputFiles/DefectiveCrystal.txt") << dc;
    fs::copy_file(root / "Library/DislocationDynamics/DD.txt", folder / "inputFiles/DD.txt");
    auto dd=read(folder / "inputFiles/DD.txt");
    parameter(dd,"coreSize","2");
    // Root minimum width is 50*sqrt(3); preserve h=8 for this fixture.
    parameter(dd,"ewaldWidthFactor","0.09237604307034013");
    parameter(dd,"climbPreconditionerCutoff","50");
    std::ofstream(folder / "inputFiles/DD.txt") << dd;
    DefectiveCrystalBase<3> base(folder.string());
    MicrostructureGenerator generator(base);
    const auto grain = base.poly.grain(1);
    const Eigen::Vector3d center = Eigen::Vector3d::Zero();
    const Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    const auto key = GlidePlaneKey<3>(center, grain->reciprocalLatticeDirection(normal));
    const auto plane = base.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(key);
    const Eigen::Vector3d burgers(0, 0, std::sqrt(8.0 / 3.0));
    const double halfWidth = rank > 0 ? 25.0 : 20.0;
    const double centerX = rank > 0 ? 40.0 : 0.0;
    std::vector<Eigen::Vector3d> polygon{{centerX - halfWidth, -10, 0},
                                         {centerX + halfWidth, -10, 0},
                                         {centerX + halfWidth, 10, 0},
                                         {centerX - halfWidth, 10, 0}};
    require(generator.insertJunctionLoop(polygon, plane, burgers, normal, center, 1, DislocationLoopIO<3>::SESSILELOOP),
            "Failed to insert climb loop");
    MicrostructureContainer<3> container(base);
    auto cd = std::make_shared<ClusterDynamics<3>>(container);
    container.push_back(cd);
    auto networkOwner=std::make_shared<DislocationNetwork<3>>(container);
    container.push_back(networkOwner);
    auto& network=*networkOwner;
    network.initializeConfiguration(generator.config(), std::ofstream(), std::ofstream());
    const auto original = network.loops().begin()->second.lock();
    const auto duplicate = network.loops().create(original->burgers(), original->glidePlane);
    std::vector<std::shared_ptr<DislocationLoopNode<3>>> nodes;
    for (const auto *link : original->linkSequence())
    {
        const auto &node = link->source;
        nodes.push_back(network.loopNodes().create(duplicate, node->networkNode, node->get_P(),
                                                   node->periodicPlanePatch(), node->periodicPlaneEdge));
    }
    network.insertLoop(duplicate, nodes);
    network.updateGeometry();

    DislocationFieldBuffers<3> buffers;
    std::size_t expected=0;
    for(const auto& entry : network.networkLinks())
    {
        const auto link=entry.second.lock();
        require(bool(link),"Expired network link");
        DislocationFieldSampling<3>::append(buffers,*link,0.2,0.1,DislocationTargetPurpose::glide);
        if(!link->hasZeroBurgers())
        {
            const auto samples=DislocationFieldSampling<3>::sample(*link,0.2);
            expected+=samples.size();
            Eigen::Matrix3d moment=Eigen::Matrix3d::Zero();
            for(std::size_t q=0;q<samples.size();++q)
            {
                const DislocationQuadraturePoint<3> legacy(link->source->sID,link->sink->sID,
                    static_cast<int>(q),static_cast<int>(samples.size()),link->sfCoeffs(),link->hermiteDofs());
                for(int axis=0;axis<3;++axis)
                {
                    require(std::abs(samples[q].position[axis]-legacy.r[axis])<1.e-12,"Position differs from existing quadrature");
                    require(std::abs(samples[q].tangent[axis]-legacy.rl[axis])<1.e-12,"Tangent differs from existing quadrature");
                }
                require(std::abs(samples[q].lineWeight-legacy.dL)<1.e-12,"Weight differs from existing quadrature");
                moment+=link->burgers()*legacy.rl.transpose()*samples[q].lineWeight;
            }
            require((moment-link->burgers()*link->chord().transpose()).norm()<1.e-10,"Integrated junction tensor");
        }
    }
    require(expected>0 && buffers.sources().size()==expected,"Missing stationary or shared sources");
    require(buffers.targets().size()==0,"Sessile loop emitted glide targets");
    require(buffers.segments().size()==network.networkLinks().size(),"Unique network-link mapping");
    auto& evaluator=network.fieldEvaluator;
    require(!evaluator.hasSnapshot(),"Evaluator must allocate lazily");
    evaluator.rebuild(0.2,0.1,DislocationTargetPurpose::glide);
    require(evaluator.buffers().sources().positions.values==buffers.sources().positions.values,"Evaluator source positions differ");
    require(evaluator.buffers().sources().weightedBurgersTensors.values==buffers.sources().weightedBurgersTensors.values,"Evaluator junction tensors differ");
    require(evaluator.buffers().targets().size()==0,"Evaluator emitted sessile glide targets");
    require(&evaluator.tree()==base.fmmContext->boxTree.get(),"Evaluator does not share FMM tree");
    bool failed=false;
    try { evaluator.rebuild(-1,0.1,DislocationTargetPurpose::glide); }
    catch(const std::invalid_argument&) { failed=true; }
    require(failed && evaluator.buffers().sources().size()==expected,"Failed rebuild altered snapshot");
    evaluator.rebuild(0.2,0.1,DislocationTargetPurpose::climb);
    require(evaluator.buffers().targets().size()>0,"Evaluator omitted climb targets");
    base.fmmContext->boxTree->setDepth(base.fmmContext->boxTree->depth()==0?1:0);
    failed=false;
    try { evaluator.leafIndex(); } catch(const std::logic_error&) { failed=true; }
    require(failed,"Evaluator accepted stale tree");
#ifdef MODELIB_TEST_EVALUATOR_KOKKOS
    evaluator.rebuild(0.2,0.1,DislocationTargetPurpose::climb);
    ScreenedDislocationOptions options;options.coreWidth=2;options.screenWidth=8;options.shortCutoff=50;
    if(rank==0)
    {
        const auto sharedSplits=network.segmentSplits();
        const auto splitBuilds=network.segmentSplitCache.buildCount();
        evaluator.prepareFMMStress(1,.3,options);
        const auto stresses=evaluator.evaluateFMMStress();
        require(stresses.stresses.values[0].size()==evaluator.buffers().targets().size(),"Missing screened stresses");
        for(const auto& component:stresses.stresses.values) {
      for(double x:component) {
        require(std::isfinite(x),"Nonfinite screened stress");
      }
    }
    for(const auto& entry:network.networkLinks()) {
      entry.second.lock()->createQuadraturePoints(true);
    }
    GalerkinClimbSolver<3> solver(network,cd.get());
        const auto snapshot=solver.screenedAssemblySnapshot();
        require(!snapshot.channels.empty()&&!snapshot.targets.empty(),"Missing network climb maps");
        if(master)
        {
            network.prepareFMMStress();
            const auto sampledStress=network.stress(Eigen::Vector3d(1,2,3),nullptr,nullptr,nullptr);
            require(sampledStress.allFinite(),"Production stress dispatch failed");
            auto probe=network.solvedQuantityState();
            for(Eigen::Index n=0;n<probe.size()/(ClusterDynamics<3>::mSize+3);++n) {
        for(int species=0;species<ClusterDynamics<3>::mSize;++species) {
          probe(n*(ClusterDynamics<3>::mSize+3)+species)=1e-10*(1+species)*(1+.1*n);
        }
      }
      network.restoreSolvedQuantityState(probe);
            cd->solvedQuantityBoundaryCache=true;
            cd->applyBoundaryConditions();
            require(cd->solvedQuantityInfluence.empty(),"Accelerated CD allocated dense H");
            require(bool(cd->solvedQuantityConcentrationEvaluator),"CD did not prepare concentration FMM");
            const auto plan=cd->solvedQuantityConcentrationEvaluator;
            const auto before=cd->clusterDynamicsFEM->mobileClusters->periodicNodalOffset();
            double error2=0,reference2=0;
            for(const auto& node:cd->clusterDynamicsFEM->mobileClusters->fe().nodes())
            {
                if(!cd->solvedQuantityTargetRows.count(node.gID)) {
          continue;
        }
        const auto direct=network.mobileConcentration(node.P0,&node,nullptr,nullptr);
                const auto actual=-before.segment(ClusterDynamics<3>::mSize*node.gID,ClusterDynamics<3>::mSize);
                error2+=(actual-direct).squaredNorm();reference2+=direct.squaredNorm();
            }
            std::cout<<"CD direct-reference relative error="<<std::sqrt(error2/reference2)<<std::endl;
            require(reference2>0&&std::sqrt(error2/reference2)<.02,"CD FMM differs from direct Helmholtz concentration");
            for(Eigen::Index n=0;n<probe.size()/(ClusterDynamics<3>::mSize+3);++n) {
        for(int species=0;species<ClusterDynamics<3>::mSize;++species) {
          probe(n*(ClusterDynamics<3>::mSize+3)+species)*=2;
        }
      }
      network.restoreSolvedQuantityState(probe);cd->applyBoundaryConditions();
            require(cd->solvedQuantityConcentrationEvaluator==plan,"Velocity change rebuilt CD FMM geometry");
            const auto after=cd->clusterDynamicsFEM->mobileClusters->periodicNodalOffset();
            require((after-2*before).norm()<1e-12*std::max(before.norm(),1e-100),"CD FMM retained stale source strengths");
            cd->solvedQuantityBoundaryCache=false;
            // Nonzero synthetic concentration drive exercises production matrix-free solve.
            for(const auto& item:network.networkLinks()) {
        for(auto& qp:item.second.lock()->quadraturePoints()){qp.cDD.setConstant(1e-8);qp.cCD.setZero();}
      }
      solver.computeClimbScalarVelocities(true);
            const auto first=solver.scalarVelocities();
            solver.computeClimbScalarVelocities(false);
            require(first.size()==solver.scalarVelocities().size(),"FMM coupled reuse lost nodes");
            for(const auto& value:solver.scalarVelocities()) {
        require(value.allFinite(),"Nonfinite FMM climb velocity");
      }
    }
        auto assembly=solver.assembleScreened(options);
        require(assembly->nodeCount()==network.networkNodes().size(),"Incorrect climb node count");
        for(std::size_t species=0;species<assembly->speciesCount();++species)
        { const auto response=assembly->apply(species,Eigen::VectorXd::Ones(assembly->nodeCount()));require(response.allFinite(),"Nonfinite network climb response"); }
        require(network.segmentSplits()==sharedSplits&&network.segmentSplitCache.buildCount()==splitBuilds,
            "Stress, CD, climb or GMRES rebuilt unchanged shared segment cuts");
        evaluator.clear();
        require(network.segmentSplits()==sharedSplits,"Evaluator clear discarded geometry split cache");
        evaluator.rebuild(.2,.1,DislocationTargetPurpose::climb);
    }
    else
    {
        bool rejected=false;try{evaluator.prepareFMMStress(1,.3,options);}catch(const std::logic_error&){rejected=true;}
        require(rejected,"Periodic screened stress must reject unsupported correction");
        GalerkinClimbSolver<3> solver(network,cd.get());rejected=false;
        try{solver.screenedAssemblySnapshot();}catch(const std::logic_error&){rejected=true;}
        require(rejected,"Periodic screened climb must reject unsupported correction");
    }
    evaluator.upload();
    evaluator.upload();
    evaluator.releaseDevice();
    evaluator.upload();
#endif
    evaluator.clear();
    require(!evaluator.hasSnapshot(),"Evaluator clear failed");
    fs::remove_all(folder);
}
int main(int argc,char** argv)
{
#ifdef MODELIB_TEST_EVALUATOR_KOKKOS
    Kokkos::ScopeGuard guard(argc,argv);
    checkEwaldInput();
#endif
    require(argc==2,"Expected repository path");
    check(argv[1],0);
    check(argv[1],1);
#ifdef MODELIB_TEST_EVALUATOR_KOKKOS
    check(argv[1],0,1);
    check(argv[1],0,2);
#endif
    std::cout<<"Real segment sampling tests passed\n";
}
