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
            ("boundary-junction-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        for(int scenario=0;scenario<6;++scenario)
        {
            const bool full=scenario==1 || scenario==2;
            {
                std::ofstream f(input/"DefectiveCrystal.txt");
                f<<"physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize="<<(full ? "0 0 0" : "0")<<";\nEwaldLengthFactor=0;\n";
                std::ofstream p(input/"polycrystal.txt");
                p<<std::setprecision(17)<<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nC2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF="<<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs="<<(full ? -1 : seeds[scenario==0 || scenario==5 ? 1 : 0])<<";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            require(base.periodicLatticeBasis.cols()==(full ? 3 : 1),"Incorrect active periodic rank");
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
            const bool yFace=scenario>=4;
            const double xFace=scenario>=2 ? -1.0 : 0.0;
            const V p0=length*(yFace ? V(-0.7,0.5,0) : V(xFace,-0.2,0));
            const V p1=length*(yFace ? V(-0.3,0.5,0) : V(xFace,0.2,0));
            std::vector<std::shared_ptr<Network::LoopType>> loops;
            for(int k=0;k<2;++k)
            {
                const auto loop=network.loops().create(burgers,plane);
                loops.push_back(loop);
                const V apex=length*(yFace ? V(-0.5,0.1+0.1*k,0) : V(-0.4-0.1*k,0,0));
                std::vector<std::shared_ptr<Network::LoopNodeType>> ring;
                for(const auto& pos : {p0,p1,apex})
                {
                    const auto node=network.networkNodes().create(pos);
                    ring.push_back(network.loopNodes().create(loop,node,pos,loop->periodicGlidePlane->getPatch(V::Zero()),
                        std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>())));
                }
                network.insertLoop(loop,ring);
            }
            network.updateBoundaryNodes();
            network.updateGeometry();
            const auto sharedLength=[&]()
            {
                double result=0;
                for(const auto& entry : network.networkLinks())
                {
                    const auto link=entry.second.lock();
                    if(link->loopLinks().size()==2)
                    {
                        result+=link->chord().norm();
                        require(std::abs(link->burgers().norm()-2*burgers.norm())<1.e-8,"Junction Burgers sum changed");
                    }
                }
                return result;
            };
            require(sharedLength()==0,"Fixture already contains shared segments");
            // Each pass can join one endpoint of initially disconnected overlaps.
            network.junctionsMaker.formJunctions(0.0);
            network.junctionsMaker.formJunctions(0.0);
            network.updateGeometry();
            const double joined=sharedLength();
            if(scenario==5)
            {
                require(joined<1.e-8,"Same-grain periodic seam received physical-boundary junction rule");
            }
            else
            {
                require(std::abs(joined-(p1-p0).norm())<1.e-7,"Physical/GB overlap did not become a shared junction");
            }
            for(const auto& loop : loops)
            {
                require(loop->isLoop(),"Junction formation broke loop closure");
                for(const auto* link : loop->linkSequence())
                {
                    const auto node=link->source;
                    require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                            "Junction formation broke periodic mapping");
                }
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: internal, periodic GB, physical-face and same-grain seam junctions\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
