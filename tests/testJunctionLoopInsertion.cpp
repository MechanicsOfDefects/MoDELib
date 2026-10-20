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
            ("junction-insertion-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
            model::MicrostructureGenerator generator(base);
            const V a=length*V(-0.5,0,0);
            const V b=length*V(-0.5,periodic ? 0.7 : 0.3,0);
            const V burgers=2.0/std::sqrt(3.0)*V::UnitY();
            const auto insert=[&](const std::vector<V>& polygon,const V& normal)
            {
                const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.get(model::GlidePlaneKey<3>(polygon.front(),grain->reciprocalLatticeDirection(normal)));
                return generator.insertJunctionLoop(polygon,plane,burgers,normal,polygon.front(),1,model::DislocationLoopIO<3>::GLISSILELOOP);
            };
            const auto snapshot=[&]()
            {
                std::ostringstream out;
                out<<std::setprecision(17);
                for(const auto& node : generator.config().nodes())
                {
                    out<<node;
                }
                for(const auto& loop : generator.config().loops())
                {
                    out<<loop;
                }
                for(const auto& node : generator.config().loopNodes())
                {
                    out<<node;
                }
                for(const auto& link : generator.config().loopLinks())
                {
                    out<<link;
                }
                for(const auto& entry : generator.uniqueNetworkNodeMap)
                {
                    out<<entry.first.transpose()<<entry.second;
                }
                return out.str();
            };
            require(insert({a,b,length*V(-0.7,0.1,0)},V::UnitZ()),"Valid first loop rejected");
            require(insert({a,b,length*V(-0.5,0.1,0.2)},V::UnitX()),"Shared junction loop rejected");
            const auto before=snapshot();
            // Internal GB, cross-grain periodic/physical back face, and a
            // nonperiodic face (Z is inactive except in the full scenario).
            require(!insert({a,length*V(0.1,0.1,0),length*V(-0.7,0.2,0)},V::UnitZ()),"Internal GB crossing accepted");
            require(!insert({a,length*V(-1.1,0.1,0),length*V(-0.7,0.2,0)},V::UnitZ()),"Back-face crossing accepted");
            if(!full)
            {
                require(!insert({a,length*V(-0.5,0.1,0.7),length*V(-0.5,0.2,0)},V::UnitX()),"Physical face crossing accepted");
            }
            require(before==snapshot(),"Rejected insertion changed configuration or sharing map");
            network.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());
            size_t shared=0;
            size_t shifted=0;
            for(const auto& entry : network.networkLinks())
            {
                const auto segment=entry.second.lock();
                if(segment->loopLinks().size()==2)
                {
                    ++shared;
                    require(segment->glidePlanes().size()==2,"Junction lost a contributing glide plane");
                    require((segment->burgers().norm()-2.0*burgers.norm())<1.e-8 && segment->burgers().norm()>burgers.norm(),"Junction Burgers sum changed");
                }
                for(const auto* plane : segment->glidePlanes())
                {
                    require(segment->source->glidePlanes().count(plane) && segment->sink->glidePlanes().count(plane),"Endpoint missing junction plane");
                }
            }
            for(const auto& node : generator.config().loopNodes())
            {
                shifted+=node.periodicShift.norm()>1.e-8;
            }
            require(shared>0,"No shared junction tested");
            require((shifted>0)==periodic,"Same-grain periodic crossing not exercised");
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: topology-based insertion and multi-plane junction sharing\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
