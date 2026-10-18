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
            ("bicrystal-open-boundary-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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

        model::DislocationDynamicsBase<3> base(folder.string());
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DDconfigIO<3> saved(folder.string());
        std::ifstream stream(root/"tests/data/bicrystalBoundary621.txt");
        require(bool(stream),"Missing bicrystal junction regression fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        using V=Eigen::Vector3d;
        const double dt=147.3858585805344;
        network.trialNodeMotion=true;
        for(const auto& entry : network.networkNodes())
        {
            const auto node=entry.second.lock();
            node->trySet_P(node->get_P()+dt*node->get_V());
        }
        const auto target=network.loops().at(162).lock();
        std::vector<V> polygon;
        for(const auto* link : target->linkSequence())
        {
            if(!link->source->isPeriodicSubdivision())
            {
                polygon.push_back(link->source->get_P());
            }
        }
        // Independent area oracle: the host grain is the slab -length<=x<=0;
        // its y/z faces are same-grain periodic. Clip only for the test oracle,
        // without choosing any walk continuation by geometric proximity.
        for(bool left : {true,false})
        {
            std::vector<V> clipped;
            for(size_t i=0; i<polygon.size(); ++i)
            {
                const V& a=polygon[i];
                const V& b=polygon[(i+1)%polygon.size()];
                const double da=left ? a.x()+length : -a.x();
                const double db=left ? b.x()+length : -b.x();
                if(da>=0)
                {
                    clipped.push_back(a);
                }
                if((da>=0)!=(db>=0))
                {
                    clipped.push_back(a+da/(da-db)*(b-a));
                }
            }
            polygon=std::move(clipped);
        }
        const auto area=[](const std::vector<V>& points)->V
        {
            V result=V::Zero();
            if(!points.empty())
            {
                for(size_t i=0; i<points.size(); ++i)
                {
                    result+=0.5*(points[i]-points.front()).cross(points[(i+1)%points.size()]-points.front());
                }
            }
            return result;
        };
        const V expected=area(polygon);
        require(expected.norm()>1,"Regression target has no resolved slipped area");
        network.updateBoundaryNodes();
        network.trialNodeMotion=false;
        const auto check=[&]()
        {
            std::vector<V> actual;
            for(const auto* link : target->linkSequence())
            {
                actual.push_back(link->source->get_P());
                require(link->sink->next.second->source==link->sink,"Broken loop connectivity");
            }
            require((area(actual)-expected).norm()<1.e-8*expected.norm(),
                    "Boundary walk lost or changed a lobe of the retained contour");
            for(const auto& entry : network.loopNodes())
            {
                const auto node=entry.second.lock();
                require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                        "Loop/network periodic mapping is inconsistent");
            }
        };
        check();
        const size_t nodes=network.loopNodes().size();
        network.updateBoundaryNodes();
        check();
        require(network.loopNodes().size()==nodes,"Repeated reconciliation changed contour topology");
        network.executeSingleGlideStepDiscreteEvents(621);
        check();
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: frame 621 open-boundary walk, area and periodic mapping\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
