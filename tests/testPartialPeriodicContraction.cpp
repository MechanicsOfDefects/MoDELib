#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
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

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("partial-contraction-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        model::StaticID<model::Polyhedron<3,3-1>>::force_count(0);

        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        const double length=200.0/std::sqrt(3.0);
        for(int scenario=0;scenario<3;++scenario)
        {
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs="<<seeds[0]<<";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            require(base.periodicLatticeBasis.cols()==1,"Incorrect active periodic rank");
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
            const V start=length*V(0.4,0,0);
            const V commonPosition=length*V(0.1,0,0);
            const auto a=network.networkNodes().create(start);
            const auto b=network.networkNodes().create(start);
            const auto common=network.networkNodes().create(commonPosition);
            std::vector<std::shared_ptr<Network::LoopType>> loops;
            const auto makeNode=[&](const auto& loop,const auto& node)
            {
                return network.loopNodes().create(loop,node,node->get_P(),loop->periodicGlidePlane->getPatch(V::Zero()),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
            };
            for(int k=0;k<3;++k)
            {
                const auto loop=network.loops().create(burgers,plane);
                loops.push_back(loop);
                const auto end=k==0 ? a : b;
                const auto apex=network.networkNodes().create(length*V(0,(k+1)*0.1,0));
                std::vector<std::shared_ptr<Network::LoopNodeType>> nodes{
                    makeNode(loop,common),makeNode(loop,end),makeNode(loop,apex)};
                if(scenario==1 && k>0)
                {
                    std::reverse(nodes.begin(),nodes.end());
                }
                network.insertLoop(loop,nodes);
            }
            network.trialNodeMotion=true;
            const V target=length*(scenario==2 ? V(0.7,0.7,0) : V(0.7,0,0));
            a->trySet_P(target);
            b->trySet_P(target);
            network.updateBoundaryNodes();
            network.trialNodeMotion=false;
            network.updateGeometry();
            if(scenario==2)
            {
                std::vector<std::shared_ptr<Network::NetworkNodeType>> physical;
                for(const auto& entry : network.networkNodes())
                {
                    const auto node=entry.second.lock();
                    for(const auto* ln : node->loopNodes())
                    {
                        if(ln->isPhysicalBoundaryNode())
                        {
                            physical.push_back(node);
                            break;
                        }
                    }
                }
                bool contracted=false;
                for(size_t i=0;i<physical.size() && !contracted;++i)
                {
                    for(size_t j=i+1;j<physical.size() && !contracted;++j)
                    {
                        if((physical[i]->get_P()-physical[j]->get_P()).norm()<1.e-8 &&
                           physical[i]->meshFaces()==physical[j]->meshFaces())
                        {
                            const auto count=physical[i]->loopNodes().size()+physical[j]->loopNodes().size();
                            require(network.contract(physical[i],physical[j]),"Physical face contraction failed");
                            require(physical[i]->loopNodes().size()+physical[j]->loopNodes().size()==count,
                                    "Physical face contraction lost loop representations");
                            contracted=true;
                        }
                    }
                }
                require(contracted,"Fixture lacks coincident physical face nodes");
                network.updateBoundaryNodes();
                // Physical clipping must preserve closure and never create a
                // periodic continuation through the inactive Y faces.
                for(const auto& entry : network.loopNodes())
                {
                    const auto node=entry.second.lock();
                    if(node->isPhysicalBoundaryNode())
                    {
                        require(!node->isPeriodicSubdivision(),"Physical corner became periodic subdivision");
                    }
                }
            }
            else
            {
                require(a->loopNodes().size()==1 && b->loopNodes().size()==2,"Fixture lost original junction sharing");
                size_t before=network.networkNodes().size();
                require(network.contract(a,b),"Endpoint contraction failed");
                require(network.networkNodes().size()+1<before,"Periodic subdivision nodes were not contracted");
                require(a->loopNodes().size()==3,"Contracted endpoint lost shared loop representations");
                size_t sharedCrossings=0;
                for(const auto& entry : network.networkNodes())
                {
                    const auto node=entry.second.lock();
                    if(node->isBoundaryNode() && node->loopNodes().size()==3)
                    {
                        ++sharedCrossings;
                    }
                }
                require(sharedCrossings>=2,"Periodic seam lost shared junction subdivisions");
                network.updateBoundaryNodes();
            }
            for(const auto& loop : loops)
            {
                require(loop->isLoop(),"Contraction broke loop closure");
                for(const auto* link : loop->linkSequence())
                {
                    const auto node=link->source;
                    require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                            "Contraction broke periodic position mapping");
                }
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: partial periodic junction contraction and physical corners\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
