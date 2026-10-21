#include <ElasticDeformation/ElasticDeformationBase.h>
#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <MicrostructureGenerator.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DislocationNetwork.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>

void require(bool value,const char* message)
{
    if(!value)
    {
        throw std::runtime_error(message);
    }
}

template<class Function>
void mustFail(Function operation,const std::string& message)
{
    bool failed=false;
    try
    {
        operation();
    }
    catch(const std::runtime_error& error)
    {
        failed=std::string(error.what()).find(message)!=std::string::npos;
    }
    require(failed,"Expected operation to fail with the specified diagnostic");
}

int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder=std::filesystem::temp_directory_path()/
            ("physics-bases-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input=folder/"inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder/"F");
        std::filesystem::copy_file(root/"Library/Materials/W.txt",input/"W.txt");
        std::filesystem::copy_file(root/"Library/Meshes/unitCube24.msh",input/"mesh.msh");
        auto configure=[&](const std::string& physics,double length,const std::string& faces)
        {
            std::ofstream(input/"DefectiveCrystal.txt")
                <<"physics="<<physics<<";\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\n"
                <<"outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
            std::ofstream(input/"polycrystal.txt")<<std::setprecision(17)
                <<"materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="
                <<length<<" 0 0 0 "<<length<<" 0 0 0 "<<length<<";\nX0=0 0 0;\nperiodicFaceIDs="<<faces<<";\n";
        };
        configure("ElasticDeformation",123.456,"-1");
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            require(base.periodicLatticeBasis.cols()==3,"Expected periodic mesh");
            require(!base.dislocationDynamicsBase && !base.clusterDynamicsBase && base.elasticDeformationBase,
                    "Active physics did not determine contexts");
            require(base.hasPhysics("ELASTICDEFORMATION"),"Case-insensitive lookup failed");
            require(&base.elasticDeformation()==base.elasticDeformationBase.get(),"Contexts are not shared");
            mustFail([&](){base.dislocationDynamics();},"not enabled");
            model::MicrostructureGenerator generator(base);
            mustFail([&](){generator.addPrismaticLoopDensity(model::PrismaticLoopDensitySpecification());},"not enabled");
            mustFail([&](){generator.insertNetworkNode(Eigen::Vector3d::Zero());},"not enabled");
            require(generator.configIO.nodes().empty() && generator.configIO.loops().empty(),"Disabled DD mutated configuration");
            model::MicrostructureContainer<3> container(base);
            mustFail([&](){model::DislocationNetwork<3> network(container);},"not enabled");
        }
        configure("DislocationDynamics",123.456,"-1");
        mustFail([&](){model::DefectiveCrystalBase<3> base(folder.string());},"not a lattice vector");
        configure("DislocationDynamics ElasticDeformation dislocationdynamics",200/std::sqrt(3.0),"-1");
        std::shared_ptr<model::DislocationDynamicsBase<3>> context;
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            require(base.activePhysics.size()==2,"Duplicate physics activation");
            require(base.dislocationDynamicsBase && base.elasticDeformationBase && !base.clusterDynamicsBase,"Incorrect mixed contexts");
            model::MicrostructureGenerator generator(base);
            model::MicrostructureContainer<3> container(base);
            require(&generator.ddBase.dislocationDynamics()==&container.ddBase.dislocationDynamics(),"Generator and physics use different contexts");
            context=base.dislocationDynamicsBase;
        }
        // A retained context owns its mesh/material even after its common base dies.
        const auto grain=context->polyOwner->grain(1);
        const auto plane=context->glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(Eigen::Vector3d::Zero(),grain->reciprocalLatticeDirection(Eigen::Vector3d::UnitZ())));
        require(bool(plane),"Context dependencies did not survive common-base destruction");
        configure("none",123.456,"");
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            require(base.activePhysics.empty() && !base.dislocationDynamicsBase && !base.clusterDynamicsBase && !base.elasticDeformationBase,"Disabled contexts were constructed");
        }
        configure("typo",123.456,"");
        mustFail([&](){model::DefectiveCrystalBase<3> base(folder.string());},"Unknown physics");
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: active physics contexts, DD validation, generator guards, and context lifetime\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
