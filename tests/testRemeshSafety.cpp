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

namespace model
{
    struct RemeshTestAccess
    {
        using Network=DislocationNetwork<3>;
        static bool eligible(Network& network,const std::pair<size_t,size_t>& key)
        {
            return network.networkRemesher.zeroChordCandidate(network.networkLinks().get(key));
        }
        static size_t contract(Network& network,const std::pair<size_t,size_t>& key)
        {
            return network.networkRemesher.contractZeroChordCandidates({key});
        }
    };
}

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("remesh-safety-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        model::StaticID<model::Polyhedron<3,3-1>>::force_count(0);
        const double length=3671.9477120460206;
        // A positive threshold distinguishes fresh area from a newly created
        // loop's initially zero cached area.
        {
            std::ifstream source(input/"DD.txt");
            std::string text((std::istreambuf_iterator<char>(source)),std::istreambuf_iterator<char>());
            const auto begin=text.find("absoluteAreaThreshold=");
            require(begin!=std::string::npos,"Missing area threshold");
            text.replace(begin,text.find(';',begin)-begin,"absoluteAreaThreshold=1");
            std::ofstream(input/"DD.txt")<<text;
        }
        for(int scenario=0;scenario<4;++scenario)
        {
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nC2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs="<<seeds[1]<<";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
            const V origin=scenario==0 ? V::Zero() : V(-0.5*length,0,0);
            const auto node=[&](const V& offset)
            {
                return network.networkNodes().create(origin+offset);
            };
            const auto insert=[&](const auto& nodes)
            {
                const auto loop=network.loops().create(burgers,plane);
                const auto patch=loop->periodicGlidePlane->getPatch(V::Zero());
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                for(const auto& n : nodes)
                {
                    std::shared_ptr<model::PeriodicPlaneEdge<3>> edge;
                    for(const auto& candidate : patch->edges())
                    {
                        if(candidate->meshIntersection->contains(n->get_P()))
                        {
                            edge=candidate;
                            break;
                        }
                    }
                    ring.push_back(network.loopNodes().create(loop,n,n->get_P(),patch,
                        std::make_pair(edge,std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
                }
                network.insertLoop(loop,ring);
                return loop;
            };
            if(scenario==0)
            {
                const auto a=node(V::Zero());
                const auto b=node(V::Zero());
                const auto first=insert(std::vector{a,b,node(V(-10,10,0)),node(V(-10,-10,0))});
                const auto second=insert(std::vector{a,b,node(V(-20,10,0)),node(V(-20,-10,0))});
                network.updateGeometry();
                require(a->isGrainBoundaryNode() && !a->isBoundaryNode(),"Expected internal GB node");
                const V probe(100,100,100);
                const Eigen::Matrix3d stress=network.stress(probe,nullptr,nullptr,nullptr);
                network.networkRemesher.contract0chordSegments();
                require(b->loopNodes().empty(),"Coincident GB nodes did not contract");
                require(a->loopNodes().size()==2 && first->isLoop() && second->isLoop(),"GB contraction lost sharing");
                network.updateGeometry();
                require((network.stress(probe,nullptr,nullptr,nullptr)-stress).norm()<1.e-10,"GB contraction changed stress");
            }
            else if(scenario==1)
            {
                const auto a=node(V::Zero());
                const auto b=node(V(20,0,0));
                const auto c=node(V(0,20,0));
                auto loop=insert(std::vector{a,b,c});
                require(loop->slippedArea()==0 && loop->currentAreaVector().norm()>100,"Expected stale initial area");
                network.networkRemesher.removeCollapsedLoops();
                require(loop->isLoop(),"Positive-area loop deleted using stale zero area");
                network.updateGeometry();
                const V oldNormal=loop->rightHandedUnitNormal();
                c->trySet_P(origin+V(0,-20,0));
                require(loop->rightHandedUnitNormal().dot(oldNormal)>0.99,"Motion overwrote normal history");
                require(loop->currentAreaVector().dot(oldNormal)<0,"Expected current contour to flip");
                network.networkRemesher.removeCollapsedLoops();
                require(loop->loopLinks().empty(),"Freshly inverted loop was not removed");
                require((loop->rightHandedUnitNormal()-oldNormal).norm()==0,"Remeshing overwrote orientation history");

                loop.reset();
                const auto shrinking=insert(std::vector{node(V(40,0,0)),node(V(60,0,0)),node(V(40,20,0))});
                network.updateGeometry();
                auto apex=shrinking->linkSequence().back()->source->networkNode;
                for(const auto* ll : shrinking->linkSequence())
                {
                    if(ll->source->get_P().y()>10)
                    {
                        apex=ll->source->networkNode;
                    }
                }
                apex->trySet_P(origin+V(40,0.01,0));
                require(shrinking->slippedArea()>100 && shrinking->currentAreaVector().norm()<1,"Expected stale large area");
                network.networkRemesher.removeCollapsedLoops();
                require(shrinking->loopLinks().empty(),"Fresh small-area loop was not removed");
            }
            else if(scenario==2)
            {
                const auto a=node(V::Zero());
                const auto b=node(V::Zero());
                const auto loop=insert(std::vector{a,b,node(V(20,20,0)),node(V(-20,20,0))});
                network.updateGeometry();
                const std::pair<size_t,size_t> key(a->sID,b->sID);
                require(model::RemeshTestAccess::eligible(network,key),"Expected initial zero-chord candidate");
                b->trySet_P(origin+V(10,0,0));
                const auto count=network.networkNodes().size();
                require(model::RemeshTestAccess::contract(network,key)==0,"Stale contraction candidate was accepted");
                require(network.networkNodes().size()==count && b->loopNodes().size()==1 && loop->isLoop(),
                        "Rejected stale contraction changed topology");
                require((b->get_P()-origin-V(10,0,0)).norm()<1.e-10,"Rejected contraction moved a node");
            }
            else
            {
                const auto a=node(V(-10,0,0));
                const auto n=node(V::Zero());
                const auto b=node(V(10,0,0));
                const auto first=insert(std::vector{a,n,b,node(V(0,10,0))});
                const auto second=insert(std::vector{b,n,node(V(-10,-10,0)),node(V(10,-10,0)),n,a,node(V(0,20,0))});
                network.updateGeometry();
                const auto before=network.loopLinks().size();
                for(const auto* ll : first->linkSequence())
                {
                    if(ll->source->networkNode==n)
                    {
                        require(!ll->source->isRemovable(25,100,0.001).first,"Distinct twin occurrences incorrectly removable");
                    }
                }
                require(network.loopLinks().size()==before,"Removal predicate changed topology");
                network.networkRemesher.remeshByRemoval();
                require(n->loopNodes().size()==3 && first->isLoop() && second->isLoop(),"Repeated visits lost junction sharing");
            }
            for(const auto& entry : network.loopNodes())
            {
                const auto ln=entry.second.lock();
                require((ln->get_P()+ln->periodicPlanePatch()->shift-ln->networkNode->get_P()).norm()<1.e-7,
                        "Remeshing broke periodic mapping");
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: internal GB contraction, fresh loop geometry, stale candidates, repeated loop visits\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
