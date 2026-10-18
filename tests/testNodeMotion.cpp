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

int main(int argc, char** argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("node-motion-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        for(int scenario=0; scenario<5; ++scenario)
        {
            const int activeAxis=(scenario==0 || scenario==3 || scenario==4) ? 0 : 1;
            {
                std::ofstream f(input / "DefectiveCrystal.txt");
                f << "physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\n"
                  << "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input / "polycrystal.txt");
                p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\n"
                  << "C2G1=1 0 0 0 1 0 0 0 1;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
                  << ";\nX0=0 0 0;\nperiodicFaceIDs=";
                if(activeAxis>=0)
                {
                    p << seeds[activeAxis];
                }
                p << ";\n";
            }
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();

            const auto first=network.loops().create(burgers,plane);
            const auto secondPlane=scenario==4 ? base.glidePlaneFactory.getFromKey(
                model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitY()))) : plane;
            const auto second=network.loops().create(2.0*burgers,secondPlane);
            const bool surface=scenario==2 || scenario==3;
            const V a=length*(scenario==3 ? V(0.4,0.5,0) : scenario==2 ? V(0.5,0,0) : V(0.4,0,0));
            const V b=length*V(0.1,0,0);
            const auto na=network.networkNodes().create(a);
            const auto nb=network.networkNodes().create(b);
            const auto makeNode=[&](const auto& loop,const auto& net,const V& position)
            {
                return network.loopNodes().create(loop,net,position,loop->periodicGlidePlane->getPatch(V::Zero()),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
            };
            network.insertLoop(first,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                makeNode(first,na,a),makeNode(first,nb,b),
                makeNode(first,network.networkNodes().create(length*V(0,0.2,0)),length*V(0,0.2,0))});
            network.insertLoop(second,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                makeNode(second,nb,b),makeNode(second,na,a),
                makeNode(second,network.networkNodes().create(length*(scenario==4 ? V(0,0,-0.2) : V(0,-0.2,0))),length*(scenario==4 ? V(0,0,-0.2) : V(0,-0.2,0)))});
            network.updateBoundaryNodes();
            const V requested=length*(scenario==2 ? V(-0.2,0.1,0) : scenario==3 ? V(0.1,-0.2,0) : V(0.3,0,0));
            na->set_V(requested,true);
            if(surface)
            {
                const int normal=scenario==2 ? 0 : 1;
                require(std::abs(na->get_V()[normal])<1.e-8,"Physical surface velocity was not projected");
            }
            network.moveNodes(1.0);
            require(!network.trialNodeMotion,"Trial mode leaked after motion");
            require(first->isLoop() && second->isLoop(),"Motion broke loop closure");
            if(scenario==1)
            {
                require(na->loopNodes().empty(),"Exited original node was not trimmed");
                bool sharedExit=false;
                for(const auto* link : first->linkSequence())
                {
                    const auto node=link->source->networkNode;
                    if(std::abs(node->get_P().x()-0.5*length)<1.e-7 && std::abs(node->get_P().y())<1.e-7)
                    {
                        sharedExit=node->loopNodes().size()==2;
                    }
                }
                require(sharedExit,"Physical exit lost junction sharing");
            }
            else
            {
                const V expected=length*(scenario==2 ? V(0.5,0.1,0) : scenario==3 ? V(0.5,0.5,0) : V(-0.3,0,0));
                require((na->get_P()-expected).norm()<1.e-7,"Incorrect committed position");
                require(na->loopNodes().size()==2,"Motion lost original junction identity");
                require((na->get_V()-(surface ? length*(scenario==2 ? V(0,0.1,0) : V(0.1,0,0)) : requested)).norm()<1.e-7,
                    "Motion discarded velocity state");
            }
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                require(base.mesh.search(node->get_P()).first,"Committed a network node outside the RVE");
            }
            for(const auto* link : first->linkSequence())
            {
                require((link->source->get_P()+link->source->periodicPlanePatch()->shift-link->source->networkNode->get_P()).norm()<1.e-7,
                    "Loop and network positions disagree");
                if(const auto segment=link->networkLink())
                {
                    if(segment->loopLinks().size()==2)
                    {
                        require(std::abs(segment->burgers().norm()-burgers.norm())<1.e-7,"Motion changed junction Burgers vector");
                    }
                }
            }
            if(scenario==3)
            {
                na->set_V(length*V(0.2,-0.1,0),true);
                require((na->get_V()-length*V(0.2,0,0)).norm()<1.e-7,"Periodic corner pinned surface sliding");
                network.moveNodes(1.0);
                require((na->get_P()-length*V(-0.3,0.5,0)).norm()<1.e-7,"Mixed-corner continuation failed");
                require(na->loopNodes().size()==2,"Mixed corner split junction ownership");
            }
            if(scenario==0)
            {
                const V oldPosition=na->get_P();
                network.trialNodeMotion=true;
                na->trySet_P(oldPosition+length*V(0.05,0,0));
                require((na->get_P()-oldPosition).norm()<1.e-8,"Trial motion prematurely committed network position");
                network.updateBoundaryNodes();
                network.trialNodeMotion=false;
                require((na->get_P()-oldPosition-length*V(0.05,0,0)).norm()<1.e-7,"Trial position not committed");
                require(na->trySet_P(na->get_P()+length*V(0.05,0,0)),"Immediate motion failed");
                require((na->get_P()-oldPosition-length*V(0.1,0,0)).norm()<1.e-7,"Standalone trySet_P did not commit");
            }
            // Repeat with zero motion to exercise retained surface anchors and seams.
            for(const auto& entry : network.networkNodes())
            {
                entry.second.lock()->set_V(V::Zero(),true);
            }
            network.moveNodes(1.0);
            require(first->isLoop() && second->isLoop(),"Repeated motion broke topology");
        }
        std::filesystem::remove_all(folder);
        std::cout << "PASS: periodic crossing, physical exit, surface sliding and mixed corner motion\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
