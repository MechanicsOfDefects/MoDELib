#include <DislocationDynamics/DislocationDynamicsBase.h>
#include <DislocationNetwork.h>
#include <DefectiveCrystal/MicrostructureContainer.h>
#include <DefectiveCrystal/DefectiveCrystalBase.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
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
            ("node-motion-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        const auto input = folder / "inputFiles";
        std::filesystem::create_directories(input);
        std::filesystem::create_directories(folder / "F");
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        const auto meshFile = root / "Library/Meshes/unitCube24.msh";
        std::filesystem::copy_file(meshFile, input / "mesh.msh");
        model::SimplicialMesh<3> mesh(meshFile.string(), Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(), {-1});
        std::vector<int> seeds;
        for(int axis = 0; axis < 3; ++axis)
        {
            for(const auto& region : mesh.regions())
            {
                for(const auto& face : region.second->faces())
                {
                    const auto& pair = face.second->periodicFacePair;
                    if(pair.second && pair.first[axis] > 0.9 &&
                       pair.first.normalized().cross(Eigen::Vector3d::Unit(axis)).norm() < 1.e-8)
                    {
                        seeds.push_back(face.first);
                    }
                }
            }
        }
        require(seeds.size() == 3, "Expected three cube periodic directions");
        model::StaticID<model::Polyhedron<3,3-1>>::force_count(0);

        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        const double length=200.0/std::sqrt(3.0);
        for(int scenario=0;scenario<5;++scenario)
        {
            {
                std::ofstream f(input / "DefectiveCrystal.txt");
                f << "physics=DislocationDynamics;\nfmmTreeDepth=0;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input / "polycrystal.txt");
                p << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF=" << std::setprecision(17) << length << " 0 0 0 " << length << " 0 0 0 " << length << ";\nX0=0 0 0;\nperiodicFaceIDs=";
                if(scenario==3)
                {
                    p << seeds[0];
                }
                p << ";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),grain->reciprocalLatticeDirection(V::UnitZ())));
            const V burgers=2.0/std::sqrt(3.0)*V::UnitX();
            const auto first=network.loops().create(burgers,plane);
            const auto second=network.loops().create(scenario==1 ? V(-burgers) : scenario==2 ? V(2*burgers) : burgers,plane);
            const auto a=network.networkNodes().create(V(0,0,0));
            const auto b=network.networkNodes().create(V(20,0,0));
            const auto c=network.networkNodes().create(V(0,20,0));
            const auto d=network.networkNodes().create(V(-20,0,0));
            const auto e=network.networkNodes().create(V(40,0,0));
            const auto makeNode=[&](const auto& loop,const auto& net)
            {
                const V offset=(scenario==3 && loop==second) ? V(length,0,0) : V::Zero();
                return network.loopNodes().create(loop,net,V(net->get_P()+offset),loop->periodicGlidePlane->getPatch(-offset),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
            };
            if(scenario==4)
            {
                network.insertLoop(first,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                    makeNode(first,c),makeNode(first,a),makeNode(first,b),makeNode(first,c),
                    makeNode(first,b),makeNode(first,e),makeNode(first,d),makeNode(first,a)});
            }
            else
            {
                network.insertLoop(first,std::vector<std::shared_ptr<Network::LoopNodeType>>{makeNode(first,c),makeNode(first,a),makeNode(first,b)});
            }
            // Opposite orientation with equal Burgers, or equal orientation
            // with opposite Burgers, both cancel the artificial CA and BC links.
            if(scenario==1)
            {
                network.insertLoop(second,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                    makeNode(second,c),makeNode(second,a),makeNode(second,d),makeNode(second,e),makeNode(second,b)});
            }
            else if(scenario!=4)
            {
                network.insertLoop(second,std::vector<std::shared_ptr<Network::LoopNodeType>>{
                    makeNode(second,a),makeNode(second,c),makeNode(second,b),makeNode(second,e),makeNode(second,d)});
            }
            const Network::LoopLinkType* link=nullptr;
            for(const auto* candidate : first->linkSequence())
            {
                if(candidate->source->networkNode==a && candidate->sink->networkNode==b)
                {
                    link=candidate;
                }
            }
            require(link!=nullptr,"Missing tested segment");
            const auto actual=link->spline();
            if(scenario==0)
            {
                // Reproduce the old loop-local choice: it bends a physically
                // straight segment toward the cancelling construction edges.
                const model::CatmullRomSplineSegment<3> oldSpline(c->get_P(),a->get_P(),b->get_P(),c->get_P());
                require(oldSpline.get_rll(0.5).norm()>1.e-3,"Regression geometry does not expose the old curvature error");
            }
            const model::CatmullRomSplineSegment<3> expected(
                scenario==2 ? c->get_P() : d->get_P(),a->get_P(),b->get_P(),scenario==2 ? c->get_P() : e->get_P());
            for(double u : {0.1,0.5,0.9})
            {
                require((actual.get_rl(u)-expected.get_rl(u)).norm()<1.e-10,"Wrong junction spline tangent");
                require((actual.get_rll(u)-expected.get_rll(u)).norm()<1.e-10,"Wrong junction spline curvature");
                if(scenario!=2)
                {
                    require(actual.get_rll(u).norm()<1.e-10,"Artificial loop closure produces line tension on a straight segment");
                }
            }
        }
        std::filesystem::remove_all(folder);
        std::cout << "PASS: cancelling junction continuation, reversed Burgers and genuine junction preservation\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
