#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <vector>
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
            ("bicrystal-corner-1915-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        const auto grain=base.poly.grain(2);
        const auto plane=base.dislocationDynamics().periodicGlidePlaneFactory.get(model::GlidePlaneKey<3>(Eigen::Vector3d(902.05206058187161,-1503.4201009697861,1202.7360807758289),grain->reciprocalLatticeDirection(Eigen::Vector3d(0.42426406871192851,-0.70710678118654757,0.56568542494923801))));
        std::vector<Eigen::Vector2d> vertices;
        std::ifstream source(root/"tests/data/bicrystalCorner1915.txt");
        require(bool(source),"Missing frame 1915 contour");
        double x,y;
        while(source>>x>>y)
        {
            vertices.emplace_back(x,y);
        }
        require(vertices.size()==6,"Incomplete frame 1915 contour");
        double referenceArea=0;
        for(size_t start=0;start<vertices.size();++start)
        {
            for(bool reverse : {false,true})
            {
                std::vector<std::pair<Eigen::Vector2d,const int* const>> polygon;
                for(size_t i=0;i<vertices.size();++i)
                {
                    polygon.emplace_back(vertices[(start+(reverse ? vertices.size()-1-i : i))%vertices.size()],nullptr);
                }
                const auto result=plane->polygonPatchIntersection(polygon,false);
                require(result.size()>=3,"Corner clipping discarded the retained area");
                double area=0;
                for(size_t i=0;i<result.size();++i)
                {
                    const Eigen::Vector2d a=std::get<0>(result[i]);
                    const Eigen::Vector2d b=std::get<0>(result[(i+1)%result.size()]);
                    area+=0.5*(a.x()*b.y()-a.y()*b.x());
                    const auto patch=plane->getPatch(std::get<1>(result[i]));
                    // Every returned node must lie inside its actual convex patch,
                    // including near the oblique periodic/physical corner.
                    for(const auto& edge : patch->edges())
                    {
                        const Eigen::Vector2d d=*edge->sink-*edge->source;
                        const Eigen::Vector2d offset=a-*edge->source;
                        require((d.x()*offset.y()-d.y()*offset.x())/d.norm()>-1.e-8,
                                "Clipping retained a node outside its patch");
                    }
                }
                if(!reverse)
                {
                    referenceArea=area;
                    require(std::abs(area)>1.0,"Corner clipping erased slipped area");
                }
                else
                {
                    require(std::abs(area+referenceArea)<1.e-6,"Corner clipping depends on polygon orientation");
                }
            }
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: frame 1915 near-coincident seam vertices in both orientations and every cyclic ordering\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<std::endl;
        return 1;
    }
}
