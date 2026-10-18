#include <DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <iostream>
#include <numbers>

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
        require(argc==2,"Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto fixture=root/"tests/data/polycrystal34";
        const auto folder=std::filesystem::temp_directory_path()/
            ("polycrystal-closed-fold-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input=folder/"inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder/"F");
        std::filesystem::create_directories(folder/"evl");
        for(const auto* name : {"polycrystal.txt","DD.txt"})
        {
            std::filesystem::copy_file(fixture/name,input/name);
        }
        std::filesystem::copy_file(root/"Library/Materials/W.txt",input/"W.txt");
        std::filesystem::copy_file(root/"Library/Meshes/poly30_periodicXYZ.msh",input/"poly30_periodicXYZ.msh");
        {
            std::ofstream out(input/"DefectiveCrystal.txt");
            out<<"physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0 0 0;\nEwaldLengthFactor=0;\n";
        }
        model::DislocationDynamicsBase<3> base(folder.string());
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DDconfigIO<3> saved(folder.string());
        std::ifstream stream(fixture/"evl_34.txt");
        require(bool(stream),"Missing frame 34 fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        using V=Eigen::Vector3d;
        const double dt=101.4090812001895;
        network.trialNodeMotion=true;
        for(const auto& entry : network.networkNodes())
        {
            const auto node=entry.second.lock();
            node->trySet_P(node->get_P()+dt*node->get_V());
        }
        const auto target=network.loops().at(579).lock();
        std::vector<V> polygon;
        for(const auto* link : target->linkSequence())
        {
            if(!link->source->isPeriodicSubdivision())
            {
                polygon.push_back(link->source->get_P());
            }
        }
        // Independent convex-grain clipping oracle for this central-patch loop.
        // Only the test uses half-space clipping; production walks port adjacency.
        for(const auto& entry : target->glidePlane->grain.region.faces())
        {
            const auto& face=entry.second;
            const V normal=face->outNormal(target->glidePlane->grain.grainID);
            std::vector<V> clipped;
            for(size_t i=0;i<polygon.size();++i)
            {
                const V& a=polygon[i];
                const V& b=polygon[(i+1)%polygon.size()];
                const double da=normal.dot(face->center()-a);
                const double db=normal.dot(face->center()-b);
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
        network.executeSingleGlideStepDiscreteEvents(34);
        check();
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: polycrystal frame 34 closed-boundary fold; frame 621 open-boundary walk, area and periodic mapping\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
