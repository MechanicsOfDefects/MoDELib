#include <FluxStreamlineTrace.h>
#include <iostream>
#include <stdexcept>
#include <cmath>
using namespace model;
void require(bool condition,const char* message)
{
    if(!condition)
    {
        throw std::runtime_error(message);
    }
}
int main(int argc,char** argv)
{
    try
    {
        require(argc==2,"Expected repository path");
        const std::string path=std::string(argv[1])+"/Library/Meshes/unitCube24.msh";
        const Eigen::Vector3d seed(0.031,0.071,0.093);
        FluxTraceOptions options;
        options.step=0.021;options.length=2.0;
        const FluxEvaluator constant=[](const Eigen::Vector3d& x,const Simplex<3,3>&)
        {
            FluxSample sample;
            sample.flux=Eigen::Vector3d::UnitX();
            return sample;
        };
        SimplicialMesh<3> freeMesh(path,Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{});
        for(const int direction : {1,-1})
        {
            const auto trace=traceFlux(freeMesh,seed,constant,options,direction);
            require(trace.size()==1 && trace[0].size()>2,"Missing physical-boundary trace");
            const auto& end=trace[0].back();
            const double boundary=direction>0?freeMesh.xMax().x():freeMesh.xMin().x();
            require(std::abs(end.position.x()-boundary)<1.e-7,"Trace did not end at physical face");
            require((end.position.tail<2>()-seed.tail<2>()).norm()<1.e-7,"Constant flux changed transverse position");
        }
        SimplicialMesh<3> periodicMesh(path,Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{-1});
        const auto periodic=traceFlux(periodicMesh,seed,constant,options,1);
        require(periodic.size()>=3,"Periodic trace did not wrap");
        double length=0.0;
        for(const auto& part : periodic)
        {
            for(size_t i=1;i<part.size();++i)
            {
                const double jump=(part[i].position-part[i-1].position).norm();
                require(jump<1.01*options.step,"Periodic wrap drawn across box");
                length+=jump;
            }
        }
        require(std::abs(length-options.length)<1.e-6,"Periodic tracing lost arclength");
        options.minimumFlux=2.0;
        require(traceFlux(freeMesh,seed,constant,options,1).front().empty(),"Minimum flux ignored");
        options.minimumFlux=0.0;
        const FluxEvaluator zero=[](const Eigen::Vector3d&,const Simplex<3,3>&){return FluxSample{};};
        require(traceFlux(freeMesh,seed,zero,options,1).front().empty(),"Zero flux generated trace");
        const FluxEvaluator rotation=[](const Eigen::Vector3d& x,const Simplex<3,3>&)
        {
            FluxSample sample;sample.flux=Eigen::Vector3d(-x.y(),x.x(),0.0);return sample;
        };
        options.length=0.3;options.step=0.002;
        const Eigen::Vector3d circleSeed(0.2,0.0,0.031);
        const auto circle=traceFlux(freeMesh,circleSeed,rotation,options,1);
        require(circle.size()==1 && circle[0].size()>20,"Missing curved trace");
        const auto last=circle[0].back().position;
        require(std::abs(last.head<2>().norm()-0.2)<2.e-5,"Midpoint integration failed curved field");
        // A periodic GB is a sink, not a seam through which flux lines wrap.
        SimplicialMesh<3> bicrystal(std::string(argv[1])+"/Library/Meshes/bicrystal_48.msh",Eigen::Matrix3d::Identity(),Eigen::Vector3d::Zero(),{-1});
        bool tested=false;
        for(const auto& pair : bicrystal.periodicSimplexPairs())
        {
            const auto* cell=pair.first->parents().begin()->second.get();
            const auto* other=pair.second.partner->parents().begin()->second.get();
            if(cell->region->regionID!=other->region->regionID)
            {
                const Eigen::Vector3d point=pair.first->center();
                const Eigen::Vector3d start=0.99*point+0.01*cell->center();
                const Eigen::Vector3d velocity=(point-start).normalized();
                const FluxEvaluator outward=[velocity](const Eigen::Vector3d&,const Simplex<3,3>&)
                {
                    FluxSample sample;sample.flux=velocity;return sample;
                };
                options.step=0.01;options.length=1.0;
                const auto trace=traceFlux(bicrystal,start,outward,options,1);
                require(trace.size()==1 && (trace[0].back().position-point).norm()<1.e-7,"Trace crossed a periodic GB sink");
                tested=true;
                break;
            }
        }
        require(tested,"Missing periodic GB fixture");
        std::cout<<"PASS: 3D curved flux, physical sinks, periodic seams, periodic GBs, forward/backward and zero flux\n";
    }
    catch(const std::exception& e)
    {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
