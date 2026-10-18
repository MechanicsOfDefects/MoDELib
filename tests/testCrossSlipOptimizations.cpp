#include <DislocationNetwork.h>
#include <DefectiveCrystal.h>
#include <CrossSlipModels.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
using Network=model::DislocationNetwork<3>;
void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
// Reference algorithms from before the cross-slip allocation optimizations.
std::pair<bool,int> legacyBase(const Network::NetworkLinkType& link,
                              const model::AthermalCrossSlipModel<Network>& model,
                              const bool enforceRemeshLength=true)
{
        if(   !link.isBoundarySegment()
           && !link.isGrainBoundarySegment()
           && link.loopLinks().size()==1
           && link.unitDirection().cross(link.burgers().normalized()).norm()<=model.sinCrossSlip
           && link.chordLength()>FLT_EPSILON
           && (!enforceRemeshLength || link.chordLength()>2.0*link.network().networkRemesher.Lmin)
           && link.grains().size()==1
           )
        {
            if((*link.loopLinks().begin())->loop->slipSystem())
            {
                const auto& grain(**link.grains().begin());
                return std::make_pair(true,grain.grainID);
            }
        }
        return std::make_pair(false,-1);
    }
std::pair<bool,std::pair<int,int>> legacySelection(const Network::NetworkLinkType& link,
                                                 const model::AthermalCrossSlipModel<Network>& model)
{

        const auto isBaseCSLink(legacyBase(link,model));

        if(isBaseCSLink.first)
        {
            // Find slip system on which screw segment has highest glide force
            const auto loop((*link.loopLinks().begin())->loop);
            const auto slipSystem(loop->slipSystem());
            if(slipSystem)
            {
                const auto& grain(**link.grains().begin());
                const auto pki=link.pkIntegral();
                const double pkGlide=std::fabs((pki-pki.dot(slipSystem->unitNormal)*slipSystem->unitNormal).norm());

                std::map<double,int> ssMap;
                for(size_t s=0;s<grain.slipSystems().size();++s)
                {
                    const auto& crossSlipSystem(grain.slipSystems().at(s));
                    if((loop->flow()-crossSlipSystem->s).squaredNorm()==0)
                    {
                        const double pkCross=std::fabs((pki-pki.dot(crossSlipSystem->unitNormal)*crossSlipSystem->unitNormal).norm());
                        ssMap.emplace(pkCross,s);
                    }
                }

                if(ssMap.size())
                {
                    const auto& crossSlipSystem(grain.slipSystems().at(ssMap.rbegin()->second));
                    if(slipSystem->unitNormal.cross(crossSlipSystem->unitNormal).squaredNorm()>FLT_EPSILON
                       && ssMap.rbegin()->first>1.1*pkGlide)
                    {// highest glide PK-force is not on current slip system
                        return std::make_pair(true,std::make_pair(grain.grainID,ssMap.rbegin()->second));
                    }
                }
            }
        }
        return std::make_pair(false,std::make_pair(-1,-1));
    }
void verify(Network& network)
{
    network.updateGeometry();
    model::AthermalCrossSlipModel<Network> model(network.ddBase.simulationParameters.traitsIO);
    for(const auto& entry : network.loops())
    {
        const auto loop=entry.second.lock();
        std::map<size_t,const Network::LoopLinkType*> ordered;
        for(const auto* link : loop->loopLinks())
        {
            ordered.emplace(link->source->sID,link);
        }
        const auto sequence=loop->linkSequence();
        require(sequence.size()==loop->loopLinks().size(),"Loop traversal count changed");
        if(!ordered.empty())
        {
            const auto* current=ordered.begin()->second;
            for(const auto* link : sequence)
            {
                require(link==current,"Loop traversal or cyclic start changed");
                current=current->next;
            }
        }
    }
    size_t eligible=0,selected=0;
    for(const auto& entry : network.networkLinks())
    {
        const auto link=entry.second.lock();
        link->createQuadraturePoints(false);
        link->updateQuadraturePoints(false);
        for(bool enforceLength : {false,true})
        {
            require(model.isBaseCrossSlipLink(*link,enforceLength)==legacyBase(*link,model,enforceLength),
                    "Cross-slip eligibility changed");
        }
        const auto legacy=legacySelection(*link,model);
        require(model.isCrossSlipLink(*link)==legacy,"Cross-slip target or tie-breaking changed");
        eligible+=legacyBase(*link,model).first;
        selected+=legacy.first;
        // Zero force exercises exact ties between every compatible plane.
        std::vector<Eigen::Vector3d> forces;
        for(auto& point : link->quadraturePoints())
        {
            forces.push_back(point.pkForce);
            point.pkForce.setZero();
        }
        require(model.isCrossSlipLink(*link)==legacySelection(*link,model),"Equal-force selection changed");
        size_t k=0;
        for(auto& point : link->quadraturePoints())
        {
            point.pkForce=forces[k++];
        }
    }
    std::cout<<"Cross-slip reference: "<<eligible<<" eligible, "<<selected<<" selected links; exact agreement\n";
}
int main(int argc, char** argv)
{
    try
    {
        require(argc == 2 || argc == 3, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        if(argc==3 && std::filesystem::is_directory(argv[2]))
        {
#ifdef _OPENMP
            omp_set_num_threads(12);
#endif
            model::DislocationDynamicsBase<3> base(argv[2]);
            model::DDconfigIO<3> saved(base.simulationParameters.traitsIO.evlFolder);
            saved.read(base.simulationParameters.runID);
            model::DefectiveCrystal<3> crystal(base);
            crystal.initializeConfiguration(saved);
            const auto network=crystal.getUniqueTypedMicrostructure<Network>();
            require(bool(network) && !network->computeDDinteractions,"Benchmark requires DD with interactions disabled");
            verify(*network);
            return 0;
        }
        const auto folder = std::filesystem::temp_directory_path() /
            ("cross-slip-optimization-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::create_directories(folder / "evl");
        const double length=3671.9477120460206;
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\n"
              << "C2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
              << ";\nX0=0 0 0;\nperiodicFaceIDs=-1;\n";
        }
        model::DislocationDynamicsBase<3> base(folder.string());
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DDconfigIO<3> saved(folder.string());
        std::ifstream stream(argc==3 ? std::filesystem::path(argv[2]) : root/"tests/data/bicrystalSharedConstraints1915.txt");
        require(bool(stream),"Missing bicrystal junction regression fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        network.computeDDinteractions=false;
#ifdef _OPENMP
        omp_set_num_threads(argc==3 ? 12 : 2);
#endif
        verify(network);
        if(argc==2)
        {
            std::shared_ptr<Network::NetworkLinkType> shared;
            for(const auto& entry : network.networkLinks())
            {
                const auto link=entry.second.lock();
                if(link->loopLinks().size()>1 && link->chordLength()>1 && !link->isBoundarySegment())
                {
                    shared=link;
                    break;
                }
            }
            require(bool(shared),"Missing shared junction for topology-change test");
            const auto memberships=shared->loopLinks().size();
            auto inserted=network.networkNodes().create(shared->get_r(0.5));
            network.expandNetworkLink(shared,inserted);
            shared.reset();
            require(inserted->loopNodes().size()==memberships,"Expansion lost junction sharing");
            verify(network);
            network.removeNetworkNode(inserted->sID);
            inserted.reset();
            verify(network);
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: exact cross-slip traversal and selection equivalence before and after shared junction topology changes\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
