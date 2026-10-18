#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <set>
#include <sstream>

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
            ("contraction-safety-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        for(int scenario=0;scenario<6;++scenario)
        {
            {
                std::ifstream original(root/"Library/DislocationDynamics/DD.txt");
                std::string dd((std::istreambuf_iterator<char>(original)),std::istreambuf_iterator<char>());
                const auto begin=dd.find("nodalVelocityConstraints=");
                require(begin!=std::string::npos,"Missing velocity constraints input");
                const auto end=dd.find(';',begin);
                dd.replace(begin,end-begin,scenario<2 ? "nodalVelocityConstraints=1 0 0" : "nodalVelocityConstraints=");
                std::ofstream(input/"DD.txt")<<dd;
            }
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs="<<seeds[0]<<";\n";
            }
            model::DislocationDynamicsBase<3> base(folder.string());
            require(!base.isPeriodicDomain,"Expected partially periodic domain");
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();

            if(scenario>=4)
            {
                const bool bent=scenario==5;
                const V pa=length*V(0.2,0,0);
                const V pb=length*V(0.2,0.1,0);
                const V pc=length*V(0.7,0,0);
                const V bend=length*V(0.5,0.2,0);
                const auto a=network.networkNodes().create(pa);
                const auto b=network.networkNodes().create(pb);
                const auto c=network.networkNodes().create(length*V(-0.3,0,0));
                const auto loop=network.loops().create(burgers,plane);
                const auto partner=network.loops().create(burgers,plane);
                const auto make=[&](const auto& owner,const auto& net,const V& position)
                {
                    return network.loopNodes().create(owner,net,position,owner->periodicGlidePlane->getPatch(V::Zero()),
                        std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
                };
                {
                    std::vector<std::shared_ptr<Network::LoopNodeType>> ring{make(loop,a,pa)};
                    if(bent)
                    {
                        ring.push_back(make(loop,network.networkNodes().create(bend),bend));
                    }
                    ring.push_back(make(loop,c,pc));
                    ring.push_back(make(loop,b,pb));
                    network.insertLoop(loop,ring);
                    const V pd=length*V(0.2,-0.2,0);
                    network.insertLoop(partner,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                        make(partner,c,pc),make(partner,a,pa),make(partner,network.networkNodes().create(pd),pd)});
                }
                network.updateBoundaryNodes();
                network.updateGeometry();
                const V partnerArea=partner->currentAreaVector();
                // Simulate an old saved configuration: seam nodes have geometry
                // but no stored distinction between bends and subdivisions.
                for(const auto& entry : network.loopNodes())
                {
                    entry.second.lock()->periodicContourVertex=false;
                }
                require(network.contract(a,b),"Periodic contour contraction was rejected");
                require((partner->currentAreaVector()-partnerArea).norm()<1.e-8,
                        "Contraction changed the shared partner's slipped area");
                require(!a->loopNodes().empty() && !c->loopNodes().empty(),
                        "Collapsed contour destroyed shared network endpoints");
                if(bent)
                {
                    require(loop->isLoop(),"Genuine periodic bend was deleted");
                    require(std::abs(loop->currentAreaVector().norm()-0.05*length*length)<1.e-8,
                            "Contraction lost slipped area at the periodic bend");
                    size_t retainedBends=0;
                    for(const auto* link : loop->linkSequence())
                    {
                        if((link->source->get_P()-bend).norm()<1.e-8)
                        {
                            require(!link->source->isPeriodicSubdivision(),"Genuine bend is still a subdivision");
                            ++retainedBends;
                        }
                    }
                    require(retainedBends>0,"Missing periodic contour vertex");
                    // Without provenance, two bulk anchors are still not a
                    // deletion certificate: retain the bent geometric contour.
                    for(const auto* link : loop->linkSequence())
                    {
                        link->source->periodicContourVertex=false;
                    }
                    network.updateBoundaryNodes();
                    require(std::abs(loop->currentAreaVector().norm()-0.05*length*length)<1.e-8,
                            "Repeated reconciliation lost the periodic bend");
                }
                else
                {
                    require(loop->linkSequence().empty(),"Retraced contour was not removed");
                    for(const auto* link : partner->linkSequence())
                    {
                        if(const auto segment=link->networkLink())
                        {
                            require(std::abs(segment->burgers().norm()-burgers.norm())<1.e-10,
                                    "Deleting a retraced loop corrupted the partner's Burgers contribution");
                        }
                    }
                }
                for(const auto& entry : network.loopNodes())
                {
                    const auto n=entry.second.lock();
                    require((n->get_P()+n->periodicPlanePatch()->shift-n->networkNode->get_P()).norm()<1.e-8,
                            "Contour contraction broke periodic mapping");
                }
                continue;
            }

            const auto snapshot=[&]()
            {
                std::ostringstream out;
                out<<std::hexfloat;
                for(const auto& entry : network.networkNodes())
                {
                    const auto n=entry.second.lock();
                    out<<n.get()<<' '<<n->sID<<' '<<n->get_P().transpose()<<' '<<n->get_V().transpose()<<'\n';
                }
                for(const auto& entry : network.loopNodes())
                {
                    const auto n=entry.second.lock();
                    out<<n.get()<<' '<<n->sID<<' '<<n->get_P().transpose()<<' '<<n->networkNode.get()<<' '
                       <<n->loop()<<' '<<n->periodicPlanePatch().get()<<' '<<n->periodicPlaneEdge.first.get()<<' '
                       <<n->prev.first<<' '<<n->next.first<<'\n';
                }
                for(const auto& entry : network.networkLinks())
                {
                    const auto link=entry.second.lock();
                    out<<link.get()<<' '<<link->source->sID<<' '<<link->sink->sID<<' '<<link->burgers().transpose()<<'\n';
                }
                for(const auto* n : network.danglingBoundaryLoopNodes)
                {
                    out<<n<<'\n';
                }
                return out.str();
            };
            const V origin=V::Zero();
            const auto a=network.networkNodes().create(origin);
            const auto b=network.networkNodes().create(scenario==2 ? V(1.e-5,0,0) : scenario<2 ? V(1,0,0) : origin);
            const auto makeNode=[&](const auto& loop,const auto& node)
            {
                return network.loopNodes().create(loop,node,node->get_P(),loop->periodicGlidePlane->getPatch(V::Zero()),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
            };
            if(scenario==3)
            {
                const auto loop=network.loops().create(burgers,plane);
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                const std::vector<V> corners{V(10,0,0),V(10,10,0),V(0,10,0),V(-10,10,0),
                                            V(-10,0,0),V(-10,-10,0),V(0,-10,0),V(10,-10,0)};
                for(size_t k=0;k<4;++k)
                {
                    ring.push_back(makeNode(loop,k%2==0 ? a : b));
                    ring.push_back(makeNode(loop,network.networkNodes().create(corners[2*k])));
                    ring.push_back(makeNode(loop,network.networkNodes().create(corners[2*k+1])));
                }
                network.insertLoop(loop,ring);
            }
            else
            {
                for(int k=0;k<2;++k)
                {
                    const auto loop=network.loops().create(burgers,plane);
                    network.insertLoop(loop,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                        makeNode(loop,k==0 ? a : b),
                        makeNode(loop,network.networkNodes().create(V(10,10+10*k,0))),
                        makeNode(loop,network.networkNodes().create(V(-10,10+10*k,0)))});
                }
            }
            network.updateGeometry();
            if(scenario<2)
            {
                const auto before=snapshot();
                require(!network.contract(scenario==0 ? a : b,scenario==0 ? b : a),"Projected target incorrectly accepted");
                require(snapshot()==before,"Rejected contraction mutated network state or junction identity");
                require(!network.contract(a,a) && snapshot()==before,"Self-contraction changed network");
            }
            else
            {
                const V probe(3,4,20);
                const auto stressBefore=network.stress(probe,nullptr,nullptr,nullptr);
                require(network.contract(a,b),"Valid contraction failed");
                require(b->loopNodes().empty(),"Removed endpoint still owns loop representations");
                for(const auto& entry : network.loopNodes())
                {
                    const auto n=entry.second.lock();
                    require((n->get_P()+n->periodicPlanePatch()->shift-n->networkNode->get_P()).norm()<1.e-10,
                            "Successful contraction did not reach exact mapped target");
                }
                network.updateGeometry();
                if(scenario==3)
                {
                    // Two A/B cuts create two additional loops. The remaining
                    // contour can still visit the surviving network node twice.
                    require(network.loops().size()==3,"Repeated visits produced incorrect loop splits");
                    require(a->loopNodes().size()==4,"Repeated-visit contraction lost junction sharing");
                    const auto stressAfter=network.stress(probe,nullptr,nullptr,nullptr);
                    require((stressAfter-stressBefore).norm()<1.e-10,"Coincident contraction changed stress/Burgers content");
                }
                for(const auto& entry : network.loops())
                {
                    require(entry.second.lock()->isLoop(),"Contraction broke loop closure");
                }
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: contraction target verification, unchanged rejection, repeated visits\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
