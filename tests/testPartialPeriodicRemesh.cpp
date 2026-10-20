#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <MicrostructureGenerator.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <set>
#include <array>

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
            ("partial-remesh-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::create_directories(folder / "config");
        model::SimplicialMesh<3> mesh((input/"mesh.msh").string(),Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{-1});
        std::array<int,3> seeds{-1,-1,-1};
        for(const auto& region : mesh.regions())
        {
            for(const auto& face : region.second->faces())
            {
                const auto& pair=face.second->periodicFacePair;
                if(pair.second)
                {
                    for(int axis=0;axis<3;++axis)
                    {
                        if(pair.first.normalized().cross(Eigen::Vector3d::Unit(axis)).norm()<1.e-8)
                        {
                            seeds[axis]=face.first;
                        }
                    }
                }
            }
        }
        require(seeds[0]>=0 && seeds[1]>=0,"Missing periodic seeds");
        model::StaticID<model::PlanarMeshFace<3>>::force_count(0);
        const double length=3671.9477120460206;
        for(int scenario=0;scenario<3;++scenario)
        {
            const bool full=scenario==1;
            const bool periodic=scenario<2;
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize="<<(full ? "0 0 0" : "0")<<";\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nC2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs=";
                if(periodic)
                {
                    p<<(full ? -1 : seeds[1]);
                }
                p<<";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
            const auto insert=[&](const auto& nodes)
            {
                const auto loop=network.loops().create(burgers,plane);
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                for(const auto& node : nodes)
                {
                    ring.push_back(network.loopNodes().create(loop,node,node->get_P(),loop->periodicGlidePlane->getPatch(V::Zero()),
                        std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
                }
                network.insertLoop(loop,ring);
                return loop;
            };
            if(periodic)
            {
                const auto a=network.networkNodes().create(length*V(-0.5,0,0));
                const auto b=network.networkNodes().create(length*V(-0.5,0.4,0));
                const auto c=network.networkNodes().create(length*V(-0.7,0.1,0));
                const auto d=network.networkNodes().create(length*V(-0.3,0.1,0));
                const auto first=insert(std::vector{a,b,c});
                const auto second=insert(std::vector{a,b,d});
                network.updateBoundaryNodes();
                network.trialNodeMotion=true;
                b->trySet_P(length*V(-0.5,0.7,0));
                network.updateBoundaryNodes();
                network.trialNodeMotion=false;
                network.updateGeometry();
                // Copy the shared arc out and back: two physical anchors plus
                // same-grain periodic subdivisions, with exactly zero area.
                std::vector<std::shared_ptr<Network::LoopNodeType>> arc;
                for(const auto* ll : first->linkSequence())
                {
                    if(ll->source->networkNode==a)
                    {
                        auto node=ll->source;
                        arc.push_back(node);
                        while(node->networkNode!=b)
                        {
                            node=network.loopNodes().get(node->next.first->sID);
                            arc.push_back(node);
                            require(arc.size()<20,"Cannot find shared periodic arc");
                        }
                        break;
                    }
                }
                require(arc.size()>2,"Missing periodic subdivisions");
                const V probe=length*V(0.3,0.2,0.1);
                const Eigen::Matrix3d stressBefore=network.stress(probe,nullptr,nullptr,nullptr);
                const auto nodeCount=network.networkNodes().size();
                auto collapsed=network.loops().create(burgers,plane);
                const auto collapsedID=collapsed->sID;
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                const auto append=[&](const auto& node)
                {
                    ring.push_back(network.loopNodes().create(collapsed,node->networkNode,node->get_P(),node->periodicPlanePatch(),node->periodicPlaneEdge));
                };
                for(const auto& node : arc)
                {
                    append(node);
                }
                for(size_t i=arc.size()-2;i>0;--i)
                {
                    append(arc[i]);
                }
                network.insertLoop(collapsed,ring);
                ring.clear();
                network.updateGeometry();
                require(collapsed->slippedArea()<FLT_EPSILON,"Fixture is not collapsed");
                collapsed.reset();
                network.networkRemesher.remove0AreaLoopAcrossBnd();
                require(!network.loops().get(collapsedID),"Collapsed periodic loop was retained");
                require(network.networkNodes().size()==nodeCount,"Cleanup lost shared network nodes");
                require(first->isLoop() && second->isLoop(),"Cleanup broke surviving loops");
                network.updateBoundaryNodes();
                network.updateGeometry();
                require((network.stress(probe,nullptr,nullptr,nullptr)-stressBefore).norm()<1.e-9*std::max(1.e-12,stressBefore.norm()),
                        "Cleanup changed shared junction stress");
                for(size_t i=1;i<arc.size();++i)
                {
                    if(arc[i]->periodicPlanePatch()==arc[i-1]->periodicPlanePatch())
                    {
                        const auto link=network.networkLinks().get(std::minmax(arc[i]->networkNode->sID,arc[i-1]->networkNode->sID));
                        require(link && link->loopLinks().size()==2 && std::abs(link->burgers().norm()-2*burgers.norm())<1.e-8,
                                "Cleanup changed junction sharing or Burgers sum");
                    }
                }
            }
            // Physical faces are real anchors even for a zero-area seed loop.
            // x=0 is an internal GB; x=-1 is a periodic GB only in the full case.
            for(const double x : {0.0,-1.0})
            {
                auto seed=network.loops().create(burgers,plane);
                const auto patch=seed->periodicGlidePlane->getPatch(V::Zero());
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                for(const double y : {-0.2,0.0,0.2})
                {
                    const V pos=length*V(x,y,0);
                    std::shared_ptr<model::PeriodicPlaneEdge<3>> edge;
                    for(const auto& candidate : patch->edges())
                    {
                        if(candidate->meshIntersection->contains(pos))
                        {
                            edge=candidate;
                            break;
                        }
                    }
                    require(bool(edge),"Missing physical edge");
                    const auto node=network.networkNodes().create(pos);
                    ring.push_back(network.loopNodes().create(seed,node,pos,patch,std::make_pair(edge,std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
                    require(ring.back()->isPhysicalBoundaryNode(),"Expected physical GB/surface anchor");
                }
                network.insertLoop(seed,ring);
                ring.clear();
                network.updateGeometry();
                require(seed->slippedArea()<FLT_EPSILON,"Expected zero-area physical seed");
                network.networkRemesher.remove0AreaLoopAcrossBnd();
                require(seed->isLoop() && seed->loopLinks().size()==3,"Cleanup removed physical boundary anchors");
            }
            for(const auto& entry : network.loopNodes())
            {
                const auto node=entry.second.lock();
                require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                        "Remeshing broke periodic mapping");
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: partial/full periodic cleanup, junction sharing, physical/GB anchors\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
