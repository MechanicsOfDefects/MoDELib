#include "ClusterPeriodicFixture.h"

// Compare endpoint flux with both the concentration derivative and direct
// integration of the differentiated kernel, including unequal endpoint rates.
void checkSegmentFlux(DislocationNetwork<3>& network,const ClusterDynamicsParameters<3>& cdp)
{
    for(const auto& entry:network.networkLinks())
    {
        const auto segment=entry.second.lock();
        if(segment->grains().size()!=1 || segment->chordLength()<1.e-6)
        {
            continue;
        }
        const Eigen::Vector3d bct=segment->burgers().cross(segment->unitDirection());
        if(std::abs(bct.dot(segment->source->climbDirection()))<1.e-8)
        {
            continue;
        }
        const auto oldA=segment->source->climbVelocityScalar;
        const auto oldB=segment->sink->climbVelocityScalar;
        segment->source->climbVelocityScalar.setConstant(0.7);
        segment->sink->climbVelocityScalar.setConstant(-0.2);
        const auto grain=(*segment->grains().begin())->grainID;
        const Eigen::Vector3d tangent=segment->unitDirection();
        const Eigen::Vector3d normal=tangent.unitOrthogonal();
        for(const auto& coordinates:std::vector<Eigen::Vector2d>{{0.4,2.0},{0.0,0.0},{0.5,0.0},{-1.0,0.01},{100.0,2.0}})
        {
            const Eigen::Vector3d x=segment->source->get_P()+coordinates.x()*segment->chord()+coordinates.y()*normal;
            const auto flux=segment->mobileClusterFlux(x,cdp);
            require(flux.allFinite(),"Nonfinite analytical segment flux");
            Eigen::Matrix<double,3,DislocationSegment<3>::mSize> gradient;
            const double step=1.e-4;
            for(int axis=0;axis<3;++axis)
            {
                const Eigen::Vector3d h=step*Eigen::Vector3d::Unit(axis);
                gradient.row(axis)=((segment->clusterConcentration(x+h,cdp)-segment->clusterConcentration(x-h,cdp))/(2*step)).transpose();
            }
            for(int k=0;k<DislocationSegment<3>::mSize;++k)
            {
                const Eigen::Vector3d fd=-cdp.D.at(grain)[k]*gradient.col(k)/cdp.omega;
                // The existing logarithmic concentration formula loses precision
                // far beyond an endpoint; use direct kernel integration there.
                if(coordinates.x()<2)
                {
                    require((fd-flux.col(k)).norm()<2.e-4*std::max(fd.norm(),flux.col(k).norm()),
                            "Analytical segment flux disagrees with concentration derivative");
                }
                double bias=1;
                if(segment->loopLinks().size()==1)
                {
                    const auto loop=(*segment->loopLinks().begin())->loop;
                    const double sign=loop->burgers().dot(loop->rightHandedUnitNormal());
                    if(std::abs(sign)>FLT_EPSILON)
                    {
                        bias=cdp.discreteDislocationBias(sign>0?0:1,k);
                    }
                }
                Eigen::Vector3d reference=Eigen::Vector3d::Zero();
                const auto& inv=cdp.invD.at(grain)[k];
                const double det=cdp.detD.at(grain)(k);
                const double core=DislocationFieldBase<3>::a2/std::cbrt(det);
                const double qA=0.7*bct.dot(segment->source->climbDirection());
                const double qB=-0.2*bct.dot(segment->sink->climbDirection());
                const int steps=4096;
                for(const auto& shift:network.ddBase.periodicShifts)
                {
                    for(int i=0;i<=steps;++i)
                    {
                        const double u=double(i)/steps;
                        const Eigen::Vector3d r=x+shift-segment->source->get_P()-u*segment->chord();
                        const double Q=r.dot(inv*r)+core;
                        const double weight=(i==0 || i==steps)?1:(i%2?4:2);
                        reference+=weight*((1-u)*qA+u*qB)*r/(Q*std::sqrt(Q));
                    }
                }
                reference*=segment->chordLength()/(3*steps*4*M_PI*std::sqrt(det)*bias*cdp.omega);
                require((reference-flux.col(k)).norm()<1.e-7*std::max(reference.norm(),flux.col(k).norm()),
                        "Analytical segment flux disagrees with integrated kernel");
            }
        }
        const auto snapshot=segment->mobileClusterFluxSource(cdp,0);
        const Eigen::Vector3d probe=segment->source->get_P()+2*normal;
        const Eigen::Vector3d saved=snapshot.evaluate(probe,network.ddBase.periodicShifts);
        require((saved-segment->mobileClusterFlux(probe,cdp).col(0)).norm()==0,"Snapshot differs from direct flux");
        segment->source->climbVelocityScalar.setZero();
        segment->sink->climbVelocityScalar.setZero();
        require(segment->mobileClusterFlux(Eigen::Vector3d::Zero(),cdp).isZero(0),"Stationary segment has nonzero flux");
        require((snapshot.evaluate(probe,network.ddBase.periodicShifts)-saved).norm()==0,"Flux snapshot changed with live source rates");
        require(segment->mobileClusterFluxSource(cdp,0).evaluate(probe,network.ddBase.periodicShifts).isZero(0),"Rebuilt snapshot retained old rates");
        segment->source->climbVelocityScalar=oldA;
        segment->sink->climbVelocityScalar=oldB;
        return;
    }
    throw std::runtime_error("No active segment for flux regression");
}

void climb(const fs::path &root, int rank)
{
    const auto folder = fixture(root, rank, 0, false);
    auto dc = read(folder / "inputFiles/DefectiveCrystal.txt");
    parameter(dc, "physics", "ClusterDynamics DislocationDynamics");
    if(rank>0)
    {
        parameter(dc,"periodicImageSize",rank==1?"1":"1 0 0");
    }
    std::ofstream(folder / "inputFiles/DefectiveCrystal.txt") << dc;
    fs::copy_file(root / "Library/DislocationDynamics/DD.txt", folder / "inputFiles/DD.txt");
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
    DislocationNetwork<3> network(container);
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
    size_t quadrature = 0;
    for (const auto &entry : network.networkLinks())
    {
        const auto link = entry.second.lock();
        link->createQuadraturePoints(true);
        for (auto &qp : link->quadraturePoints())
        {
            qp.cCD.setZero();
            qp.cDD.setZero();
            qp.cDD(0) = 1.e-4;
            ++quadrature;
        }
    }
    require(quadrature > 0, "Climb regression has no active quadrature points");
    require(bool(network.climbSolver), "Climb solver not enabled");
    network.climbSolver->computeClimbScalarVelocities(true);
    checkSegmentFlux(network,cd->cdp);
    const auto velocity = network.climbSolver->getNodeVelocities();
    require(velocity.allFinite() && velocity.norm() > 0, "Climb drive produced no finite motion");
    size_t relations = 0;
    forEachPeriodicVelocityRelation(
        network,
        [&](size_t n, size_t a, size_t b, double wa, double wb)
        {
            ++relations;
            require((velocity.segment<3>(3 * n) - wa * velocity.segment<3>(3 * a) - wb * velocity.segment<3>(3 * b))
                            .norm() < 1.e-9 * velocity.norm(),
                    "Actual climb solver violates periodic vector interpolation");
        });
    require(rank == 0 || relations > 0, "Periodic climb test has no crossings");
    size_t tangentialMarkers = 0;
    size_t localIndex = 0;
    for (const auto &entry : network.networkNodes())
    {
        const auto node = entry.second.lock();
        const auto direction = node->climbDirection();
        const auto value = velocity.segment<3>(3 * localIndex);
        const auto &species = network.climbSolver->scalarVelocities()[localIndex];
        const double normalSpeed = -(species * cd->cdp.msVector / cd->cdp.msVector.abs()).sum();
        require(std::abs(value.dot(direction) - normalSpeed) < 1.e-9 * velocity.norm(),
                "Climb source strength differs from normal motion");
        if ((value - direction * value.dot(direction)).norm() > 1.e-8 * velocity.norm())
        {
            ++tangentialMarkers;
        }
        ++localIndex;
    }
    require(rank == 0 || tangentialMarkers > 0, "Regression does not exercise tangential subdivision motion");
    network.climbSolver->computeClimbScalarVelocities(false);
    require((network.climbSolver->getNodeVelocities() - velocity).norm() < 1.e-10 * velocity.norm(),
            "Cached climb solve changed velocities");
    size_t index = 0;
    for (const auto &entry : network.networkNodes())
    {
        entry.second.lock()->set_V(velocity.segment<3>(3 * index++), true);
    }
    network.moveNodes(1.e-4 / velocity.norm());
    for (const auto &entry : network.networkLinks())
    {
        require(entry.second.lock()->loopLinks().size() == 2, "Climb motion broke shared junction ownership");
    }
    fs::remove_all(folder);
}
// Controlled thermal drive: vacancy loops emit and interstitial loops absorb.
// Test source coupling independently of the outer FEM boundary iteration.
std::map<int,double> mixedAnnealing(const fs::path& root,int selection)
{
    const auto folder=fixture(root,0,0,false);
    auto dcText=read(folder / "inputFiles/DefectiveCrystal.txt");
    parameter(dcText,"physics","ClusterDynamics DislocationDynamics");
    std::ofstream(folder / "inputFiles/DefectiveCrystal.txt")<<dcText;
    fs::copy_file(root / "Library/DislocationDynamics/DD.txt",folder / "inputFiles/DD.txt");
    DefectiveCrystalBase<3> base(folder.string());
    MicrostructureGenerator generator(base);
    const auto grain=base.poly.grain(1);
    const Eigen::Vector3d normal=Eigen::Vector3d::UnitZ();
    const auto key=GlidePlaneKey<3>(Eigen::Vector3d::Zero(),grain->reciprocalLatticeDirection(normal));
    const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(key);
    for(int sign : {-1,1})
    {
        if(selection!=0 && selection!=sign)
        {
            continue;
        }
        const Eigen::Vector3d center(20.0*sign,0,0);
        std::vector<Eigen::Vector3d> polygon;
        for(int i=0;i<12;++i)
        {
            const double angle=2.0*M_PI*i/12;
            polygon.emplace_back(center+Eigen::Vector3d(8*std::cos(angle),8*std::sin(angle),0));
        }
        require(generator.insertJunctionLoop(polygon,plane,sign*std::sqrt(8.0/3.0)*normal,
                                              normal,center,1,DislocationLoopIO<3>::SESSILELOOP),
                "Failed to insert annealing regression loop");
    }
    MicrostructureContainer<3> container(base);
    auto cd=std::make_shared<ClusterDynamics<3>>(container);
    container.push_back(cd);
    DislocationNetwork<3> network(container);
    network.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());
    const double equilibrium=cd->cdp.equilibriumMobileConcentration(0.0)(0);
    for(const auto& entry:network.networkLinks())
    {
        const auto link=entry.second.lock();
        const auto loop=(*link->loopLinks().begin())->loop;
        const double sign=loop->burgers().dot(loop->rightHandedUnitNormal())>0?1.0:-1.0;
        link->createQuadraturePoints(true);
        for(auto& qp:link->quadraturePoints())
        {
            qp.cDD.setZero();
            qp.cCD.setZero();
            qp.cCD(0)=equilibrium;
            qp.cDD(0)=equilibrium*std::exp(0.1*sign);
        }
    }
    network.climbSolver->computeClimbScalarVelocities(true);
    checkSegmentFlux(network,cd->cdp);
    const Eigen::VectorXd velocity=network.climbSolver->getNodeVelocities();
    std::map<int,double> means;
    std::map<int,int> counts;
    for(const auto& entry:network.networkNodes())
    {
        const auto node=entry.second.lock();
        const auto loop=(*node->loopNodes().begin())->loop();
        const int sign=loop->burgers().dot(loop->rightHandedUnitNormal())>0?1:-1;
        const Eigen::Vector3d center(20.0*sign,0,0);
        means[sign]+=(node->get_P()-center).normalized().dot(velocity.segment<3>(3*node->gID()));
        ++counts[sign];
    }
    for(auto& entry:means)
    {
        entry.second/=counts.at(entry.first);
        require(entry.second<0,"Thermal annealing loop must shrink");
    }
    fs::remove_all(folder);
    return means;
}
int main(int argc, char **argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        for (int rank : {0, 1, 3})
        {
            climb(argv[1], rank);
        }
        const auto vacancy=mixedAnnealing(argv[1],1);
        const auto interstitial=mixedAnnealing(argv[1],-1);
        const auto mixed=mixedAnnealing(argv[1],0);
        require(mixed.at(1)<vacancy.at(1),"Interstitial absorption must accelerate vacancy-loop shrinkage");
        require(mixed.at(-1)<interstitial.at(-1),"Vacancy emission must accelerate interstitial-loop shrinkage");
        std::cout << "PASS: mixed annealing and coupled CD/climb solver, periodic crossings, motion and junction sharing\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
