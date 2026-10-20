#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
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
        require(argc == 2 || argc == 3, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const bool frame283=argc==3 && std::string(argv[2])=="283";
        const auto folder = std::filesystem::temp_directory_path() /
            ("bicrystal-junction-collapse-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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

        model::DefectiveCrystalBase<3> base(folder.string());
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DefectiveCrystalConfigIO<3> saved(folder.string());
        std::ifstream stream(root/(frame283 ? "tests/data/bicrystalJunction283.txt" : "tests/data/bicrystalJunction596.txt"));
        require(bool(stream),"Missing bicrystal junction regression fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        // Saved velocities and the original timestep reproduce the junction
        // contraction that reduces a triangular loop to two distinct points.
        network.moveNodes(frame283 ? 272.4170214605611 : 236.5277449432792);
        network.executeSingleGlideStepDiscreteEvents(frame283 ? 283 : 596);
        network.updateBoundaryNodes();
        for(const auto& entry : network.loops())
        {
            const auto loop=entry.second.lock();
            const auto links=loop->linkSequence();
            size_t anchors=0;
            for(const auto* link : links)
            {
                if(!link->source->isPeriodicSubdivision())
                {
                    ++anchors;
                }
            }
            require(anchors>=3,"Collapsed periodic contour survived junction reconciliation");
            require(links.size()>=3,"Collapsed loop survived junction reconciliation");
            for(const auto* link : links)
            {
                require(link->sink->next.second->source==link->sink,"Loop connectivity is broken");
            }
        }
        for(const auto& entry : network.loopNodes())
        {
            const auto node=entry.second.lock();
            require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                    "Junction reconciliation broke the loop/network mapping");
        }
        // Reconciliation is idempotent, including the memberships and Burgers
        // sums of shared network links after removing a collapsed loop.
        std::map<std::pair<size_t,size_t>,Eigen::Vector3d> burgers;
        for(const auto& entry : network.networkLinks())
        {
            burgers.emplace(entry.first,entry.second.lock()->burgers());
        }
        const auto loopCount=network.loops().size();
        network.updateBoundaryNodes();
        require(network.loops().size()==loopCount,"Repeated update changed surviving loops");
        require(network.networkLinks().size()==burgers.size(),"Repeated update changed junction links");
        for(const auto& entry : network.networkLinks())
        {
            require((entry.second.lock()->burgers()-burgers.at(entry.first)).norm()<1.e-12,
                    "Repeated update changed a shared Burgers contribution");
        }
        // A near-coincident edge must be treated consistently with the walk's
        // edge tolerance; a resolvable thin triangle must remain a contour.
        const auto grain=base.poly.grain(1);
        const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.get(model::GlidePlaneKey<3>(Eigen::Vector3d::Zero(),grain->reciprocalLatticeDirection(Eigen::Vector3d::UnitZ())));
        for(double separation : {1.e-12,1.e-5})
        {
            std::vector<std::pair<Eigen::Vector2d,const int* const>> polygon;
            const Eigen::Vector2d center=plane->referencePlane->localPosition(Eigen::Vector3d(-0.5*length,0,0));
            polygon.emplace_back(center,nullptr);
            polygon.emplace_back(center+Eigen::Vector2d(separation,0),nullptr);
            polygon.emplace_back(center+Eigen::Vector2d(0,10),nullptr);
            if(separation<1.e-8)
            {
                require(plane->polygonPatchIntersection(polygon,false).empty(),
                        "Near-coincident triangle returned a two-point contour");
            }
            else
            {
                require(plane->polygonPatchIntersection(polygon,false).size()==3,
                        "Resolvable thin triangle was discarded");
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: junction collapse, shared Burgers contributions and periodic mapping\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
