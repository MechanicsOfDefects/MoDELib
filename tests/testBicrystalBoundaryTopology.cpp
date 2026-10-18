#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <map>
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
            ("bicrystal-boundary-topology-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        for(int frame : {430,433})
        {
            model::DislocationDynamicsBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            model::DislocationNetwork<3> network(container);
            model::DDconfigIO<3> saved(folder.string());
            std::ifstream stream(root/"tests/data"/("bicrystalBoundary"+std::to_string(frame)+".txt"));
            require(bool(stream),"Missing bicrystal boundary regression fixture");
            saved.readTxtStream(stream);
            network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
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
                require(after.size()==7 && after.size()==before.size(),"Boundary update split or lost a loop");
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
            network.updateBoundaryNodes();
            check(initial,1.e-10);
            network.updateBoundaryNodes();
            check(initial,1.e-10);
            // The recorded velocities reproduce the 430 branch jump and the
            // 433 non-closing walk; no field recomputation or tutorial files.
            network.trialNodeMotion=true;
            for(const auto& entry : network.networkNodes())
            {
                const auto node=entry.second.lock();
                node->trySet_P(node->get_P()+16.95870980823875*node->get_V());
            }
            const auto trial=areas();
            network.updateBoundaryNodes();
            network.trialNodeMotion=false;
            check(trial,1.e-3);
            const auto clipped=areas();
            network.updateBoundaryNodes();
            check(clipped,1.e-10);
            const auto plasticArea=[&]()
            {
                Eigen::Matrix3d result=Eigen::Matrix3d::Zero();
                for(const auto& entry : areas())
                {
                    result+=network.loops().at(entry.first).lock()->burgers()*entry.second.transpose();
                }
                return result;
            };
            const Eigen::Matrix3d plasticBefore=plasticArea();
            network.executeSingleGlideStepDiscreteEvents(frame);
            require((plasticArea()-plasticBefore).norm()<1.e-3*std::max(1.0,plasticBefore.norm()),
                    "Junction formation changed Burgers-weighted slipped area");
            // GB junctions can legitimately split loops. Check the retained
            // oriented area and mapping instead of requiring unchanged IDs.
            V areaBefore=V::Zero();
            V areaAfter=V::Zero();
            for(const auto& entry : clipped)
            {
                areaBefore+=entry.second;
            }
            for(const auto& entry : areas())
            {
                areaAfter+=entry.second;
            }
            require((areaAfter-areaBefore).norm()<1.e-3*std::max(1.0,areaBefore.norm()),
                    "Junction formation changed retained oriented area");
            for(const auto& entry : network.loops())
            {
                require(entry.second.lock()->isLoop(),"Junction formation broke loop closure");
            }
            for(const auto& entry : network.loopNodes())
            {
                const auto node=entry.second.lock();
                require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                        "Junction formation broke periodic mapping");
            }
            std::cout<<"PASS: bicrystal frame "<<frame<<" topology, mapping, area and remeshing\n";
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
