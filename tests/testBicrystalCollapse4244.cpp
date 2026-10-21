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
#include <array>
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
            ("bicrystal-collapse-4244-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        std::ifstream stream(root/"tests/data/bicrystalCollapse4244.txt");
        require(bool(stream),"Missing frame 4244 fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        require(network.loops().size()==6,"Incomplete junction neighborhood");
        require(network.loops().at(5860).lock()->currentAreaVector().norm()>1.0,"Initial loop has no slipped area");
        {
            const auto loop=network.loops().at(5860).lock();
            const auto plane=loop->periodicGlidePlane;
            const Eigen::Vector3d originalApex=network.loopNodes().get(6036331)->get_P();
            const std::array<Eigen::Vector3d,3> trial{
                Eigen::Vector3d(3672.8571855755372,480.38042416509757,-1119.2670014528765),
                Eigen::Vector3d(3671.947712046021,458.39269828388035,-1146.0695536572607),
                Eigen::Vector3d(3671.947712046021,489.69414114074192,-1106.9427500861839)};
            for(bool outside : {false,true})
            {
                for(size_t start=0;start<3;++start)
                {
                    for(bool reverse : {false,true})
                    {
                        std::vector<std::pair<Eigen::Vector2d,const int* const>> polygon;
                        for(size_t k=0;k<3;++k)
                        {
                            const size_t i=(start+(reverse ? 3-k : k))%3;
                            polygon.emplace_back(plane->referencePlane->localPosition(i==0 && !outside ? originalApex : trial[i]),nullptr);
                        }
                        const auto contour=plane->polygonPatchIntersection(polygon,false);
                        if(outside)
                        {
                            require(contour.empty(),"Clipped line returned a degenerate contour");
                        }
                        else
                        {
                            require(contour.size()>=3,"Valid nearby slipped area was discarded");
                            double area=0;
                            for(size_t k=0;k<contour.size();++k)
                            {
                                const auto a=std::get<0>(contour[k])-std::get<0>(contour.front());
                                const auto b=std::get<0>(contour[(k+1)%contour.size()])-std::get<0>(contour.front());
                                area+=a.x()*b.y()-a.y()*b.x();
                            }
                            require(std::abs(area)>1.0,"Nearby contour lost its nonzero area");
                        }
                    }
                }
            }
        }
        network.moveNodes(108.4240043273397);
        require(network.loops().count(5860)==0,"Collapsed retained line survived as a loop");
        // Three coincident cap loops collapse together. Their base nodes
        // belong only to those caps and must disappear; the adjoining three
        // nonzero-area loops must keep their shared network segments.
        const std::set<size_t> survivingLoops{5852,5892,5893};
        require(network.loops().size()==survivingLoops.size(),"Unexpected surviving loop count");
        for(const auto& entry : network.loops())
        {
            require(survivingLoops.count(entry.first) && entry.second.lock()->currentAreaVector().norm()>1.0,
                    "Lost a nonzero-area neighbor");
        }
        require(!network.networkLinks().empty(),"Neighboring junction vanished");
        for(const auto& entry : network.networkLinks())
        {
            const auto link=entry.second.lock();
            std::set<size_t> memberships;
            Eigen::Vector3d sum=Eigen::Vector3d::Zero();
            for(const auto* ll : link->loopLinks())
            {
                memberships.insert(ll->loop->sID);
                sum+=(ll->source->networkNode==link->source ? 1.0 : -1.0)*ll->flow().cartesian();
            }
            require(memberships==survivingLoops,"Clipping damaged the three-loop junction sharing");
            require((sum-link->burgers()).norm()<1.e-12,"Shared Burgers sum is inconsistent");
        }
        const auto loopCount=network.loops().size();
        const auto nodeCount=network.networkNodes().size();
        std::map<std::pair<size_t,size_t>,Eigen::Vector3d> burgers;
        for(const auto& entry : network.networkLinks())
        {
            burgers.emplace(entry.first,entry.second.lock()->burgers());
        }
        network.updateBoundaryNodes();
        require(network.loops().size()==loopCount && network.networkNodes().size()==nodeCount &&
                network.networkLinks().size()==burgers.size(),"Boundary update is not idempotent");
        for(const auto& entry : network.networkLinks())
        {
            require((entry.second.lock()->burgers().array()==burgers.at(entry.first).array()).all(),"Shared Burgers vector changed");
        }
        for(const auto& entry : network.loopNodes())
        {
            const auto node=entry.second.lock();
            require((node->get_P()+node->periodicPlanePatch()->shift-node->networkNode->get_P()).norm()<1.e-7,
                    "Periodic loop/network mapping changed");
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: frame 4244 collapsed clipping result and shared boundary junctions\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
