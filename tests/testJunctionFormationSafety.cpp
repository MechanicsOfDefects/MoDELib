#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <set>

void require(bool condition, const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}

namespace model
{
    struct JunctionFormationTestAccess
    {
        using Network=DislocationNetwork<3>;
        using Formation=DislocationJunctionFormation<Network>;
        using Queue=std::deque<Formation::IntersectionTypeContainerType>;
        static Queue find(Network& network)
        {
            Queue result(1);
            network.junctionsMaker.findIntersections(result,1);
            return result;
        }
        static size_t contract(Network& network,const Queue& queue)
        {
            return network.junctionsMaker.contractJunctions(queue);
        }
        static size_t refresh(Network& network,const Queue& queue)
        {
            size_t count=0;
            for(const auto& item : queue.front())
            {
                const auto a=network.networkLinks().get(std::get<0>(item));
                const auto b=network.networkLinks().get(std::get<1>(item));
                if(a && b)
                {
                    Queue refreshed(1);
                    network.junctionsMaker.findPairIntersections(refreshed,a.get(),b.get());
                    count+=refreshed.front().size();
                }
            }
            return count;
        }
        static void glissile(Network& network)
        {
            network.junctionsMaker.glissileJunctions(0.0);
        }
    };
}

int main(int argc, char** argv)
{
    try
    {
        require(argc >= 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("junction-safety-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        const auto meshFile = root / "Library/Meshes/unitCube24.msh";
        std::filesystem::copy_file(meshFile, input / "mesh.msh");
        model::SimplicialMesh<3> mesh(meshFile.string(), Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(), {-1});
        std::vector<int> seeds;
        for(int axis = 0; axis < 3; ++axis)
        {
            for(const auto& region : mesh.regions())
            {
                for(const auto& face : region.second->faces())
                {
                    const auto& pair = face.second->periodicFacePair;
                    if(pair.second && pair.first[axis] > 0.9 &&
                       pair.first.normalized().cross(Eigen::Vector3d::Unit(axis)).norm() < 1.e-8)
                    {
                        seeds.push_back(face.first);
                    }
                }
            }
        }
        require(seeds.size() == 3, "Expected three cube periodic directions");
        model::StaticID<model::PlanarMeshFace<3>>::force_count(0);

        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        const double length=200.0/std::sqrt(3.0);
        for(int scenario=0;scenario<5;++scenario)
        {
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs="<<seeds[scenario==2 ? 0 : 1]<<";\n";
            }
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            using Access=model::JunctionFormationTestAccess;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto planeAt=[&](const V& point,const V& normal)
            {
                return base.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(point,grain->reciprocalLatticeDirection(normal)));
            };
            const auto plane=planeAt(V::Zero(),V::UnitZ());
            const auto insert=[&](const V& burgers,const auto& gp,const auto& nodes)
            {
                const auto loop=network.loops().create(burgers,gp);
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                for(const auto& node : nodes)
                {
                    ring.push_back(network.loopNodes().create(loop,node,node->get_P(),loop->periodicGlidePlane->getPatch(V::Zero()),
                        std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
                }
                network.insertLoop(loop,ring);
                return loop;
            };
            if(scenario==0 || scenario>=3)
            {
                const double delta=8.0/std::sqrt(3.0);
                const V burgers=2.0/std::sqrt(3.0)*V::UnitY();
                const auto a=network.networkNodes().create(V(0,-20,0));
                const auto b=network.networkNodes().create(V(0,20,0));
                const auto c=network.networkNodes().create(V(-20,0,0));
                const auto d=network.networkNodes().create(V(delta,-20,0));
                const auto e=network.networkNodes().create(V(delta,20,0));
                const auto f=network.networkNodes().create(V(delta+10,0,-10));
                const auto first=insert(burgers,plane,std::vector{a,b,c});
                const auto second=insert(-burgers,planeAt(V(delta,0,0),V(1,0,1)),std::vector{d,e,f});
                network.updateGeometry();
                const auto keyA=std::make_pair(a->sID,b->sID);
                const auto keyB=std::make_pair(d->sID,e->sID);
                std::set<std::pair<size_t,size_t>> independentKeys;
                if(scenario>=3)
                {
                    // A second, independent collision must remain eligible in the
                    // same batch. A shared loop at the first collision must survive.
                    insert(2*burgers,plane,std::vector{a,b,network.networkNodes().create(V(-30,0,0))});
                    const V shift=scenario==3 ? V(0,0,40.0/std::sqrt(3.0)) : V(0,40.0/std::sqrt(3.0),0);
                    const auto aa=network.networkNodes().create(a->get_P()+shift);
                    const auto bb=network.networkNodes().create(b->get_P()+shift);
                    const auto cc=network.networkNodes().create(c->get_P()+shift);
                    const auto dd=network.networkNodes().create(d->get_P()+shift);
                    const auto ee=network.networkNodes().create(e->get_P()+shift);
                    const auto ff=network.networkNodes().create(f->get_P()+shift);
                    if(scenario==3)
                    {
                        insert(burgers,planeAt(shift,V::UnitZ()),std::vector{aa,bb,cc});
                    }
                    else
                    {
                        // Sharing even a distant node couples the two collisions.
                        insert(burgers,plane,std::vector{aa,bb,cc,c});
                    }
                    insert(-burgers,planeAt(V(delta,0,0)+shift,V(1,0,1)),std::vector{dd,ee,ff});
                    independentKeys.insert(std::minmax(aa->sID,bb->sID));
                    independentKeys.insert(std::minmax(dd->sID,ee->sID));
                    network.updateGeometry();
                }
                auto queue=Access::find(network);
                Access::Queue selected(1);
                for(const auto& item : queue.front())
                {
                    if((std::get<0>(item)==keyA && std::get<1>(item)==keyB) ||
                       (std::get<0>(item)==keyB && std::get<1>(item)==keyA) ||
                       (independentKeys.count(std::get<0>(item)) && independentKeys.count(std::get<1>(item))))
                    {
                        selected.front().push_back(item);
                    }
                }
                require(!selected.front().empty(),"X-separated attractive segments were missed by sweep");
                if(scenario>=3)
                {
                    require(selected.front().size()>=2,"Expected an event for each independent collision");
                    // Queue the remaining collisions on the same contours after
                    // the two independent pairs. These must wait for reconciliation.
                    for(const auto& event : queue.front())
                    {
                        selected.front().push_back(event);
                    }
                    require(Access::contract(network,selected)==(scenario==3 ? 2 : 1),
                            "Batch must contract independent pairs and defer their dependent events");
                    for(const auto& entry : network.loops())
                    {
                        require(entry.second.lock()->isLoop(),"Batch broke shared loop closure");
                    }
                    bool shared=false;
                    for(const auto& entry : network.networkLinks())
                    {
                        shared=shared || entry.second.lock()->loopLinks().size()>1;
                    }
                    require(shared,"Batch lost junction sharing");
                    require(Access::contract(network,Access::find(network))>0,
                            "Deferred junction was not available in the next iteration");
                }
                else
                {
                // Keep the link IDs alive but invalidate their stored collision.
                network.trialNodeMotion=true;
                for(const auto& node : {d,e,f})
                {
                    node->trySet_P(node->get_P()+V(20,0,-20));
                }
                network.updateBoundaryNodes();
                network.trialNodeMotion=false;
                network.updateGeometry();
                const auto nodeCount=network.networkNodes().size();
                const auto linkCount=network.networkLinks().size();
                require(Access::contract(network,selected)==0,"Stale collision was contracted");
                require(network.networkNodes().size()==nodeCount && network.networkLinks().size()==linkCount,
                        "Stale rejection changed topology");
                }
            }
            else
            {
                // Two sessile loops share a segment whose sum is a slip Burgers vector.
                V burgers=V::Zero();
                for(const auto& slip : grain->slipSystems())
                {
                    if(std::abs(slip.second->n.cartesian().normalized().x())<1.e-8)
                    {
                        burgers=slip.second->s.cartesian();
                        break;
                    }
                }
                require(burgers.norm()>0,"Missing X-directed glide plane");
                const auto a=network.networkNodes().create(V(0,0,0));
                const auto b=network.networkNodes().create(V(0.4*length,0,0));
                const auto c=network.networkNodes().create(V(0.2*length,0.2*length,0));
                const auto d=network.networkNodes().create(V(0.2*length,-0.2*length,0));
                const auto first=insert(2*burgers,plane,std::vector{a,b,c});
                const auto second=insert(-burgers,plane,std::vector{a,b,d});
                network.updateBoundaryNodes();
                network.updateGeometry();
                network.trialNodeMotion=true;
                b->trySet_P(V(0.7*length,0,0));
                network.updateBoundaryNodes();
                network.trialNodeMotion=false;
                network.updateGeometry();
                std::shared_ptr<Network::NetworkNodeType> endpoint;
                for(const auto& entry : network.networkLinks())
                {
                    const auto link=entry.second.lock();
                    if(link->loopLinks().size()==2 && (link->source==a || link->sink==a))
                    {
                        endpoint=link->source==a ? link->sink : link->source;
                    }
                }
                require(bool(endpoint) && endpoint->isBoundaryNode(),"Missing boundary junction endpoint");
                // A third loop visits the endpoint but does not use the shared segment.
                const auto e=network.networkNodes().create(endpoint->get_P()+V(-10,5,0));
                const auto f=network.networkNodes().create(endpoint->get_P()+V(-10,-5,0));
                const auto unrelated=insert(burgers,plane,std::vector{endpoint,e,f});
                network.updateBoundaryNodes();
                network.updateGeometry();
                const auto area=[&]()
                {
                    V result=V::Zero();
                    for(const auto* ll : unrelated->linkSequence())
                    {
                        result+=0.5*ll->source->get_P().cross(ll->sink->get_P());
                    }
                    return result;
                };
                const V before=area();
                const auto count=network.loops().size();
                const V probe(7,11,19);
                const Eigen::Matrix3d stressBefore=network.stress(probe,nullptr,nullptr,nullptr);
                Access::glissile(network);
                network.updateGeometry();
                const Eigen::Matrix3d stressAfter=network.stress(probe,nullptr,nullptr,nullptr);
                require((stressAfter-stressBefore).norm()<1.e-9*std::max(1.e-12,stressBefore.norm()),
                        "Glissile seed changed the original stress field");
                const auto shared=network.networkLinks().get(std::minmax(a->sID,endpoint->sID));
                require(shared && shared->loopLinks().size()==3 && shared->hasZeroBurgers(),
                        "Glissile seed lost the original shared segment or its Burgers cancellation");
                require(network.loops().size()>count,"Expected glissile junction was not created");
                require((area()-before).norm()<1.e-8,"Unrelated loop at shared endpoint was altered");
                require(first->isLoop() && second->isLoop() && unrelated->isLoop(),"Junction sharing broke loop closure");
            }
            for(const auto& entry : network.loopNodes())
            {
                const auto ln=entry.second.lock();
                require((ln->get_P()+ln->periodicPlanePatch()->shift-ln->networkNode->get_P()).norm()<1.e-7,
                        "Junction formation broke periodic mapping");
            }
        }
        if(argc>2 && std::string(argv[2])=="benchmark")
        {
            std::filesystem::copy_file(root/"Library/Meshes/bicrystal_48.msh",input/"mesh.msh",std::filesystem::copy_options::overwrite_existing);
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                const double scale=3671.9477120460206;
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nC2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF="<<scale<<" 0 0 0 "<<scale<<" 0 0 0 "<<scale<<";\nX0=0 0 0;\nperiodicFaceIDs=-1;\n";
            }
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            model::DislocationNetwork<3> network(container);
            model::DDconfigIO<3> saved(folder.string());
            std::ifstream stream(root/"tests/data/bicrystalBoundary433.txt");
            require(bool(stream),"Missing collision benchmark fixture");
            saved.readTxtStream(stream);
            network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
            using Access=model::JunctionFormationTestAccess;
            auto queue=Access::find(network);
            const auto searchStart=std::chrono::steady_clock::now();
            for(int i=0;i<20;++i)
            {
                queue=Access::find(network);
            }
            const double searchSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-searchStart).count()/20;
            const auto refreshStart=std::chrono::steady_clock::now();
            size_t refreshed=0;
            for(int i=0;i<200;++i)
            {
                refreshed+=Access::refresh(network,queue);
            }
            const double refreshSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-refreshStart).count()/200;
            require(!queue.front().empty() && refreshed>0,"Benchmark has no queued collisions");
            std::cout<<"BENCHMARK: "<<queue.front().size()<<" queued events; search="<<searchSeconds
                     <<" s; refresh="<<refreshSeconds<<" s; ratio="<<refreshSeconds/searchSeconds<<'\n';
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: X-separated candidates, stale collisions, physical endpoints and shared periodic junctions\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
