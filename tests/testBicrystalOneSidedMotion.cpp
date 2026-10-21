#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <OneSidedMotion.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <iterator>
void require(bool b,const char* message)
{
    if(!b)
    {
        throw std::runtime_error(message);
    }
}
int main(int argc,char** argv)
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

        using V=Eigen::Vector3d;
        using Network=model::DislocationNetwork<3>;
        // A corner with nonorthogonal physical faces must satisfy both inequalities.
        const std::vector<V> outward{V::UnitX(),V(1,1,0).normalized()};
        const std::vector<V> equality{V::UnitZ()};
        const V corner=model::oneSidedMotion<3>(V(-1,2,3),equality,outward);
        require((corner-V(0,2,0)).norm()<1.e-10,"Incorrect corner active set");
        const V blocked=model::oneSidedMotion<3>(V(-2,-1,0),equality,outward);
        require(blocked.norm()<1.e-10,"Corner allowed inward motion");
        // Opposite normals from two host grains permit only common sliding.
        require((model::oneSidedMotion<3>(V(1,2,0),equality,{V::UnitX(),-V::UnitX()})-V(0,2,0)).norm()<1.e-10,
                "Opposite host-grain normals were merged");
        for(int side : {0,1})
        {
            for(bool mixedCorner : {false,true})
            {
                model::DefectiveCrystalBase<3> base(folder.string());
                model::MicrostructureContainer<3> container(base);
                Network network(container);
                const auto grain=base.poly.grain(1);
                const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
                const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
                const auto first=network.loops().create(burgers,plane);
                const auto second=network.loops().create(2.0*burgers,plane);
                const double sign=side==0?1.0:-1.0;
                const double wall=side==0?0.0:-length;
                const double top=mixedCorner?0.5*length:0.2*length;
                const std::vector<V> positions{{wall,-0.2*length,0},{wall,top,0},
                    {wall-sign*0.3*length,top,0},{wall-sign*0.3*length,-0.2*length,0}};
                {
                    std::vector<std::shared_ptr<Network::LoopNodeType>> a,b;
                    for(const auto& p : positions)
                    {
                        const auto node=network.networkNodes().create(p);
                        const auto edges=std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>());
                        a.push_back(network.loopNodes().create(first,node,p,first->periodicGlidePlane->getPatch(V::Zero()),edges));
                        b.push_back(network.loopNodes().create(second,node,p,second->periodicGlidePlane->getPatch(V::Zero()),edges));
                    }
                    network.insertLoop(first,a);
                    network.insertLoop(second,b);
                }
                network.updateBoundaryNodes();
                network.updateGeometry();
                size_t boundary=0;
                std::vector<std::shared_ptr<Network::NetworkNodeType>> exiting;
                for(const auto& entry : network.networkNodes())
                {
                    const auto node=entry.second.lock();
                    node->set_V(V::Zero(),true);
                    if(std::abs(node->get_P().x()-wall)<1.e-7)
                    {
                        ++boundary;
                        require(node->isGrainBoundaryNode(),"Missing GB node");
                        const V old=node->get_P();
                        node->set_V(V(-sign,0,0),true);
                        require(node->get_V().norm()<1.e-10,"Inward GB velocity was allowed");
                        node->trySet_P(old+V(-sign,0,0));
                        require((node->get_P()-old).norm()<1.e-8,"Position commit allowed inward detachment");
                        const V requested(length*V(sign*0.1,0.05,0));
                        node->set_V(requested,true);
                        require((node->get_V()-requested).norm()<1.e-8,"Outward trial velocity was reduced");
                        exiting.push_back(node);
                    }
                }
                require(boundary>=2,"Missing physical boundary chain");
                // Trial movement must affect all loop representations, not network P.
                network.trialNodeMotion=true;
                for(const auto& node : exiting)
                {
                    const V old=node->get_P();
                    node->trySet_P(old+node->get_V());
                    require((node->get_P()-old).norm()<1.e-10,"Trial motion committed outside grain");
                    for(const auto* ln : node->loopNodes())
                    {
                        require((ln->get_P()+ln->periodicPlanePatch()->shift-old-node->get_V()).norm()<1.e-7,
                                "Junction representations did not share outward trial motion");
                    }
                }
                network.updateBoundaryNodes();
                network.trialNodeMotion=false;
                network.updateGeometry();
                for(const auto& node : exiting)
                {
                    require(node->loopNodes().empty(),"Outside trial anchor was retained after trimming");
                }
                require(first->isLoop() && second->isLoop(),"Trimming broke a loop");
                size_t gbNodes=0;
                for(const auto& entry : network.networkNodes())
                {
                    const auto node=entry.second.lock();
                    if(node->loopNodes().empty())
                    {
                        continue;
                    }
                    require(base.mesh.search(node->get_P()).first,"Network node outside RVE");
                    require(node->get_P().x()<=1.e-7 && node->get_P().x()>=-length-1.e-7,"Network node left host grain");
                    require(node->loopNodes().size()==2,"Retrimming split a shared junction");
                    gbNodes+=node->isGrainBoundaryNode();
                }
                require(gbNodes>=2,"Lost GB chain after trimming");
                for(const auto& entry : network.networkLinks())
                {
                    const auto link=entry.second.lock();
                    require(link->loopLinks().size()==2,"Replacement segment lost loop sharing");
                    require(std::abs(link->burgers().norm()-3.0*burgers.norm())<1.e-8,"Changed junction Burgers vector");
                    if(link->isGrainBoundarySegment())
                    {
                        require(link->meshLocation()==2,"GB segment lost physical classification");
                    }
                }
                std::cout<<"PASS side "<<side<<" mixed corner "<<mixedCorner<<std::endl;
            }
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
