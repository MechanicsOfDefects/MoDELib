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
            ("periodic-fields-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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

        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        const double length=200.0/std::sqrt(3.0);
        for(int activeAxis : {-1,2})
        {
            const bool periodic=activeAxis==0;
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
            const std::array<V,4> positions{length*V(-.8,-.8,0),length*V(.8,-.8,0),length*V(.8,.8,0),length*V(-.8,.8,0)};
            std::vector<std::shared_ptr<Network::NetworkNodeType>> owners;
            for(const auto& p : positions)
            {
                owners.push_back(network.networkNodes().create(0.5*p));
            }
            std::vector<std::shared_ptr<Network::LoopType>> loops;
            for(size_t family=0;family<3;++family)
            {
                const auto loop=network.loops().create(burgers,plane);
                loops.push_back(loop);
                std::vector<std::shared_ptr<Network::LoopNodeType>> nodes;
                for(size_t k=0;k<positions.size();++k)
                {
                    const size_t index=family==1 ? positions.size()-1-k : k;
                    const auto net=family==2 ? network.networkNodes().create(0.5*positions[index]) : owners[index];
                    nodes.push_back(network.loopNodes().create(loop,net,positions[index],loop->periodicGlidePlane->getPatch(V::Zero()),
                        std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
                }
                network.insertLoop(loop,nodes);
            }
            network.updateBoundaryNodes();
            std::set<size_t> sharedIDs;
            for(const auto* link : loops[0]->linkSequence())
            {
                const auto segment=link->networkLink();
                require(segment && segment->loopLinks().size()==2,"Enclosure lost shared contour ownership");
                require(segment->burgers().norm()<1.e-12,"Opposite enclosed loops no longer cancel");
                sharedIDs.insert(link->source->networkNode->sID);
            }
            require(sharedIDs.size()==4,"Incorrect enclosed face contour");
            for(const auto* link : loops[2]->linkSequence())
            {
                require(!sharedIDs.count(link->source->networkNode->sID),"Unrelated enclosing loop acquired junction sharing");
                require(link->networkLink()->loopLinks().size()==1,"Unrelated loop was merged");
            }
            // The reconciliation stage must retain cancelling lobes too, not
            // only the geometry-only walk exercised in testPeriodicPolygonWalk.
            const auto balanced=network.loops().create(burgers,plane);
            std::vector<std::shared_ptr<Network::LoopNodeType>> balancedNodes;
            const std::array<V,4> balancedPositions{length*V(-.2,-.2,0),length*V(.2,.2,0),length*V(-.2,.2,0),length*V(.2,-.2,0)};
            for(const V& position : balancedPositions)
            {
                balancedNodes.push_back(network.loopNodes().create(balanced,network.networkNodes().create(position),position,
                    balanced->periodicGlidePlane->getPatch(V::Zero()),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
            }
            network.insertLoop(balanced,balancedNodes);
            network.updateBoundaryNodes();
            require(balanced->isLoop() && balanced->linkSequence().size()==4,"Cancelling lobes were deleted during reconciliation");
            const auto nodeCount=network.networkNodes().size();
            network.updateBoundaryNodes();
            require(network.networkNodes().size()==nodeCount,"Enclosed boundary update is not idempotent");
            for(const auto* link : loops[0]->linkSequence())
            {
                require(sharedIDs.count(link->source->networkNode->sID),"Repeated update changed enclosed junction identity");
            }
        }
        std::cout<<"PASS: enclosure preserves shared ownership and separates unrelated loops\n";
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
    return 0;
}
