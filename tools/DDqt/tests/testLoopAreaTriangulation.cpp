#include <LoopAreaTriangulation.h>
#include <Polygon2D.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <iterator>

void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
using V=Eigen::Vector2d;
double area(const std::vector<V>& polygon)
{
    double result=0.0;
    for(const auto& t : model::triangulateLoopArea(polygon))
    {
        const V a=t[1]-t[0];
        const V b=t[2]-t[0];
        result+=0.5*(a.x()*b.y()-a.y()*b.x());
        require(model::Polygon2D::windingNumber(V((t[0]+t[1]+t[2])/3.0),polygon)!=0,
                "Triangle filled a zero-winding region");
    }
    return result;
}
void example(std::vector<V> polygon,double expected)
{
    for(int direction=0;direction<2;++direction)
    {
        for(size_t i=0;i<polygon.size();++i)
        {
            require(std::abs(area(polygon)-expected)<1.e-10,"Incorrect filled area or start-dependent tessellation");
            std::rotate(polygon.begin(),polygon.begin()+1,polygon.end());
        }
        std::reverse(polygon.begin(),polygon.end());
    }
}
int main(int argc,char** argv)
{
    try
    {
        example({{0,0},{2,0},{2,2},{0,2}},4);
        example({{0,0},{2,2},{0,2},{2,0}},2); // bow tie: zero signed area, two visible lobes
        example({{0,0},{2,0},{2,1},{1,1},{1,2},{0,2}},3);
        example({{0,0},{1,0},{1,1},{0,1},{0,0},{2,0},{3,0},{3,1},{2,1},{2,0}},2);
        example({{0,0},{3,0},{3,3},{0,3},{0,0},{1,1},{1,2},{2,2},{2,1},{1,1}},8);
        example({{0,0},{1,0},{1,1},{0,1},{0,0},{1,0},{1,1},{0,1}},1); // winding two, fill once
        example({{0,0},{1,0},{1,1},{0,1},{0,0},{0,1},{1,1},{1,0}},0);
        example({{0,0},{0,0},{1,0},{1,1},{0,1},{0,0}},1);
        example({{0,0},{1,0},{2,0}},0);
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

        for(int frame : {430,433,500})
        {
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            model::DislocationNetwork<3> network(container);
            model::DefectiveCrystalConfigIO<3> saved(folder.string());
            std::ifstream stream(root/"tests/data"/("bicrystalBoundary"+std::to_string(frame)+".txt"));
            require(bool(stream),"Missing bicrystal fixture");
            saved.readTxtStream(stream);
            network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
            size_t patches=0;
            for(const auto& item : network.loops())
            {
                const auto loop=item.second.lock();
                Eigen::Vector3d fullArea=Eigen::Vector3d::Zero();
                for(const auto* link : loop->linkSequence())
                {
                    fullArea+=0.5*link->source->get_P().cross(link->sink->get_P());
                }
                // Compare against the actual loop, not merely against the clipped
                // input to the renderer: omitted physical corners corrupt that input.
                require((loop->patches().orientedArea()-fullArea).norm()<1.e-8*std::max(1.0,fullArea.norm()),
                        "Patch assembly discarded part of the original loop area");
                for(const auto& patch : loop->patches().localPatches())
                {
                    const auto& polygon=patch.second;
                    const double original=area(polygon);
                    auto reversed=polygon;
                    std::reverse(reversed.begin(),reversed.end());
                    require(std::abs(area(reversed)-original)<1.e-8*std::max(1.0,original),
                            "Bicrystal fill depends on contour orientation");
                    const auto triangles=model::triangulateLoopArea(polygon);
                    if(polygon.empty())
                    {
                        continue;
                    }
                    V low=polygon.front(),high=low;
                    for(const auto& p : polygon)
                    {
                        low=low.cwiseMin(p);
                        high=high.cwiseMax(p);
                    }
                    // Independently compare raster coverage against the contour winding.
                    for(int i=0;i<19;++i)
                    {
                        for(int j=0;j<19;++j)
                        {
                            const V p=low+(high-low).cwiseProduct(V((i+0.371)/19,(j+0.613)/19));
                            int hits=0;
                            for(const auto& t : triangles)
                            {
                                hits+=model::Polygon2D::windingNumber(p,std::vector<V>{t[0],t[1],t[2]})!=0;
                            }
                            require(hits==(model::Polygon2D::windingNumber(p,polygon)!=0),
                                    "Bicrystal tessellation has missing or overlapping fill");
                        }
                    }
                    ++patches;
                }
            }
            require(patches>0,"No bicrystal patches tested");
            std::cout<<"Checked frame "<<frame<<": "<<patches<<" patches\n";
        }
        std::filesystem::remove_all(folder);
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return 1;
    }
}
