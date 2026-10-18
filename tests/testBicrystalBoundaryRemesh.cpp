#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <map>
#include <set>
#include <iterator>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("bicrystal-boundary-remesh-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        // Match the tutorial remeshing length, without depending on its inputs.
        {
            const auto dd=input/"DD.txt";
            std::ifstream source(dd);
            std::string contents((std::istreambuf_iterator<char>(source)),std::istreambuf_iterator<char>());
            const auto location=contents.find("Lmax=100;");
            require(location!=std::string::npos,"Expected library remeshing default");
            contents.replace(location,std::string("Lmax=100;").size(),"Lmax=150;");
            std::ofstream output(dd);
            output<<contents;
        }
        for(int frame : {430,433,500})
        {
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            model::DislocationNetwork<3> network(container);
            model::DDconfigIO<3> saved(folder.string());
            std::ifstream stream(root/"tests/data"/("bicrystalBoundary"+std::to_string(frame)+".txt"));
            require(bool(stream),"Missing bicrystal boundary regression fixture");
            saved.readTxtStream(stream);
            network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
            // Overlay one real bicrystal loop on the same network nodes to
            // exercise a shared boundary chain (including nonzero junction Burgers).
            {
                std::shared_ptr<model::DislocationNetwork<3>::LoopType> source;
                size_t mostPhysical=0;
                for(const auto& entry : network.loops())
                {
                    const auto candidate=entry.second.lock();
                    size_t physical=0;
                    for(const auto* link : candidate->linkSequence())
                    {
                        physical+=link->source->isPhysicalBoundaryNode();
                    }
                    if(physical>mostPhysical)
                    {
                        mostPhysical=physical;
                        source=candidate;
                    }
                }
                require(bool(source),"Missing physical-boundary loop");
                const auto overlay=network.loops().create(source->burgers(),source->glidePlane);
                std::vector<std::shared_ptr<model::DislocationNetwork<3>::LoopNodeType>> nodes;
                for(const auto* link : source->linkSequence())
                {
                    const auto node=link->source;
                    nodes.push_back(network.loopNodes().create(overlay,node->networkNode,node->get_P(),
                        node->periodicPlanePatch(),node->periodicPlaneEdge));
                }
                network.insertLoop(overlay,nodes);
            }
            network.updateGeometry();
            using V=Eigen::Vector3d;
            const auto areas=[&]()
            {
                std::map<size_t,V> result;
                for(const auto& entry : network.loops())
                {
                    V area=V::Zero();
                    for(const auto* link : entry.second.lock()->linkSequence())
                    {
                        area+=0.5*link->source->get_P().cross(link->sink->get_P());
                    }
                    result.emplace(entry.first,area);
                }
                return result;
            };
            const auto check=[&](const std::map<size_t,V>& before,double tolerance)
            {
                const auto after=areas();
                require(after.size()==8 && after.size()==before.size(),"Boundary update split or lost a loop");
                for(const auto& entry : before)
                {
                    require((after.at(entry.first)-entry.second).norm()<tolerance*std::max(1.0,entry.second.norm()),
                            "Boundary walk jumped to another contour branch");
                }
                for(const auto& entry : network.loopNodes())
                {
                    const auto node=entry.second.lock();
                    require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                            "Loop/network periodic mapping is inconsistent");
                }
            };
            const auto initial=areas();
            std::set<size_t> protectedNodes;
            size_t sharedBefore=0;
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                bool protect=node->neighbors().size()!=2;
                bool physical=false;
                for(const auto* ln : node->loopNodes())
                {
                    protect|=ln->isPeriodicSubdivision() || ln->periodicPlaneEdge.second!=nullptr;
                    physical|=ln->isPhysicalBoundaryNode();
                }
                if(protect)
                {
                    protectedNodes.insert(entry.first);
                }
                sharedBefore+=physical && node->loopNodes().size()>1;
            }
            const auto sharedLengths=[&]()
            {
                std::map<std::set<size_t>,double> result;
                for(const auto& entry : network.networkLinks())
                {
                    const auto segment=entry.second.lock();
                    std::set<size_t> loops;
                    for(const auto* link : segment->loopLinks())
                    {
                        loops.insert(link->loop->sID);
                    }
                    result[loops]+=segment->chordLength();
                }
                return result;
            };
            const auto stress=[&]()
            {
                Eigen::Matrix3d result=Eigen::Matrix3d::Zero();
                const V point(0.31*length,0.23*length,0.17*length);
                for(const auto& entry : network.networkLinks())
                {
                    const auto segment=entry.second.lock();
                    result+=segment->straight.stress(point);
                }
                return result;
            };
            const auto lengthsBefore=sharedLengths();
            const auto stressBefore=stress();
            const auto countBefore=network.networkNodes().size();
            const auto removed=network.networkRemesher.simplifyPhysicalBoundaryChains();
            require(removed>0,"No redundant boundary nodes removed");
            require(network.networkNodes().size()+removed==countBefore,"Network nodes leaked after removal");
            require(network.networkRemesher.simplifyPhysicalBoundaryChains()==0,"Simplification did not reach a fixed point");
            check(initial,1.e-10);
            for(const auto id : protectedNodes)
            {
                require(network.networkNodes().find(id)!=network.networkNodes().end(),"Removed a corner, crossing or branching junction");
            }
            size_t sharedAfter=0;
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                bool physical=false;
                for(const auto* ln : node->loopNodes())
                {
                    physical|=ln->isPhysicalBoundaryNode();
                }
                sharedAfter+=physical && node->loopNodes().size()>1;
            }
            require(sharedAfter<sharedBefore,"Did not simplify any shared physical-boundary chains");
            const auto lengthsAfter=sharedLengths();
            require(lengthsBefore.size()==lengthsAfter.size(),"Changed participating loop sets");
            for(const auto& entry : lengthsBefore)
            {
                require(std::abs(lengthsAfter.at(entry.first)-entry.second)<1.e-8*std::max(1.0,entry.second),
                        "Changed junction sharing or segment lengths");
            }
            network.updateGeometry();
            require((stress()-stressBefore).norm()<1.e-8*std::max(1.e-12,stressBefore.norm()),
                    "Simplification changed the network stress field");
            network.updateBoundaryNodes();
            check(initial,1.e-10);
            require(network.networkRemesher.simplifyPhysicalBoundaryChains()==0,"Boundary update recreated redundant nodes");
            std::cout<<"Frame "<<frame<<": removed "<<removed<<" network nodes, shared physical nodes "
                     <<sharedBefore<<" -> "<<sharedAfter<<std::endl;
        }
        std::filesystem::remove_all(folder);
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
