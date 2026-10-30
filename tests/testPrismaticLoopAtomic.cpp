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
        require(argc == 2 || argc == 3, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto folder = std::filesystem::temp_directory_path() /
            ("prismatic-atomic-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");

        std::filesystem::copy_file(root / "Library/Meshes/bicrystal_48.msh", input / "mesh.msh");
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        std::filesystem::create_directories(folder / "config");
        const double length=3671.9477120460206;
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
            std::ofstream p(input / "polycrystal.txt");
            p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\n"
              << "C2G2=0.8 0 0.6 0 1 0 -0.6 0 0.8;\nF=" << length << " 0 0 0 " << length << " 0 0 0 " << length
              << ";\nX0=0 0 0;\nperiodicFaceIDs=-1;\n";
        }
        model::DefectiveCrystalBase<3> base(folder.string());
        using V=Eigen::Vector3d;
        model::MicrostructureGenerator generator(base);
        const auto snapshot=[&]()
        {
            std::ostringstream out;
            out<<std::setprecision(17);
            for(const auto& x : generator.config().nodes())
            {
                out<<x<<'\n';
            }
            for(const auto& x : generator.config().loops())
            {
                out<<x<<'\n';
            }
            for(const auto& x : generator.config().loopNodes())
            {
                out<<x<<'\n';
            }
            for(const auto& x : generator.config().loopLinks())
            {
                out<<x<<'\n';
            }
            for(const auto& x : generator.uniqueNetworkNodeMap)
            {
                out<<x.first.transpose()<<' '<<x.second<<'\n';
            }
            return out.str();
        };
        const auto insert=[&](const V& center,const double steps)
        {
            model::PrismaticLoopIndividualSpecification spec;
            spec.slipSystemIDs={0};
            spec.loopRadii={1.e-7};
            spec.glideSteps={steps};
            spec.loopCenters.resize(1,3);
            spec.loopCenters.row(0)=center;
            generator.addPrismaticLoopIndividual(spec);
        };
        const auto reject=[&](const V& center,const double steps)
        {
            const auto before=snapshot();
            bool rejected=false;
            try
            {
                insert(center,steps);
            }
            catch(const std::runtime_error&)
            {
                rejected=true;
            }
            require(rejected,"Prism crossing into another grain was accepted");
            require(snapshot()==before,"Failed prism changed loops, nodes, links or the node map");
        };
        // A failed first attempt must leave every output container empty.
        reject(V(-0.02*length,0,0),50);
        require(generator.config().loops().empty() && generator.config().nodes().empty() &&
                generator.config().loopNodes().empty() && generator.config().loopLinks().empty() &&
                generator.uniqueNetworkNodeMap.empty(),"Failed first prism left residual topology");
        insert(V(-0.5*length,0,0),50);
        require(generator.config().loops().size()==7,"Contained prism is incomplete");
        // Reject both the internal GB and the cross-grain periodic back face.
        reject(V(-0.02*length,0,0),50);
        reject(V(-0.98*length,0,0),50);
        // The base fits, but the translated cap and sides cross a GB.
        reject(V(-0.5*length,0,0),10000);
        reject(V(-0.5*length,0,0),-10000);
        const size_t oldLoops=generator.config().loops().size();
        // Cross two same-grain periodic directions at a box corner.
        insert(V(-0.5*length,0.49*length,0.49*length),50);
        require(generator.config().loops().size()==oldLoops+7,"Same-grain periodic prism is incomplete");
        size_t shifted=0;
        std::map<std::pair<size_t,size_t>,V> burgers;
        std::set<std::pair<size_t,size_t>> baseEdges;
        for(const auto& node : generator.config().loopNodes())
        {
            shifted+=node.periodicShift.norm()>1.e-8;
        }
        require(shifted>0,"Periodic acceptance test never crossed a seam");
        for(const auto& link : generator.config().loopLinks())
        {
            if(link.hasNetworkLink)
            {
                const auto& source=generator.config().loopNodes().at(link.sourceID);
                const auto& sink=generator.config().loopNodes().at(link.sinkID);
                const auto& loop=generator.config().loops().at(link.loopID);
                const auto key=std::make_pair(std::min(source.networkNodeID,sink.networkNodeID),
                                              std::max(source.networkNodeID,sink.networkNodeID));
                auto entry=burgers.emplace(key,V::Zero()).first;
                entry->second+=(source.networkNodeID<sink.networkNodeID? 1.0 : -1.0)*loop.B;
                if(loop.loopType==model::DislocationLoopIO<3>::SESSILELOOP)
                {
                    baseEdges.insert(key);
                }
            }
        }
        require(!baseEdges.empty(),"No sessile base edges tested");
        for(const auto& key : baseEdges)
        {
            require(burgers.at(key).norm()<1.e-10,"Prism base edge lacks Burgers cancellation");
        }
        reject(V(-0.02*length,0,0),50);
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: atomic prism rejection, same-grain periodic crossings and base cancellation\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
