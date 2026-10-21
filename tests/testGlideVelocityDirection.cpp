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
        std::filesystem::copy_file(root / "Library/DislocationDynamics/DD.txt", input / "DD.txt");
        const double length=200.0/std::sqrt(3.0);
        {
            {
                std::ofstream f(input / "DefectiveCrystal.txt");
                f << "physics=DislocationDynamics;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=0;\nEwaldLengthFactor=0;\n";
                std::ofstream p(input / "polycrystal.txt");
                p << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF=" << std::setprecision(17) << length << " 0 0 0 " << length << " 0 0 0 " << length << ";\nX0=0 0 0;\nperiodicFaceIDs=";
                p << ";\n";
            }
            model::DefectiveCrystalBase<3> base(folder.string());
            model::MicrostructureContainer<3> container(base);
            using Network=model::DislocationNetwork<3>;
            using V=Eigen::Vector3d;
            Network network(container);
            const auto grain=base.poly.grain(1);

            const auto ss=grain->slipSystems().at(0);
            const V burgers=ss->s.cartesian();
            const V normal=ss->unitNormal;
            const V tangent=normal.cross(burgers).normalized();
            const V glide=normal.cross(tangent).normalized();
            const auto plane=base.dislocationDynamics().glidePlaneFactory.getFromKey(model::GlidePlaneKey<3>(V::Zero(),ss->n));
            const auto loop=network.loops().create(burgers,plane);
            const auto makeNode=[&](const V& p)
            {
                return network.loopNodes().create(loop,network.networkNodes().create(p),p,loop->periodicGlidePlane->getPatch(V::Zero()),
                    std::make_pair(std::shared_ptr<model::PeriodicPlaneEdge<3>>(),std::shared_ptr<model::PeriodicPlaneEdge<3>>()));
            };
            const auto a=makeNode(-10*tangent);
            const auto b=makeNode(10*tangent);
            network.insertLoop(loop,std::vector<std::shared_ptr<Network::LoopNodeType>>{a,b,makeNode(10*glide)});
            network.updateGeometry();
            const auto segment=a->next.second->networkLink();
            require(bool(segment->slipSystem()),"Test segment has no mobility slip system");
            using Q=model::DislocationQuadraturePoint<3>;
            Q point;
            const Eigen::Matrix3d stress=0.1*(burgers.normalized()*normal.transpose()+normal*burgers.normalized().transpose());
            const auto velocity=[&](const V& force)
            {
                return Q::getGlideVelocity(*segment,V::Zero(),force,stress,tangent,1.0,-1);
            };
            require(velocity(tangent).norm()<1.e-12,"Purely tangential force drives quadrature motion");
            require(velocity(normal).norm()<1.e-12,"Climb force drives glide motion");
            const V reference=velocity(glide);
            require(reference.norm()>1.e-12,"Mobility test is vacuous: no glide response");
            const V mixed=velocity(glide+100*tangent+20*normal);
            require((mixed-reference).norm()<1.e-10,"Tangential force changes glide velocity");
            require(std::abs(mixed.dot(tangent))<1.e-12,"Velocity is not orthogonal to the segment");
            require(std::abs(mixed.dot(normal))<1.e-12,"Velocity leaves the glide plane");
            require((velocity(-glide)+reference).norm()<1.e-10,"Glide force sign is lost");
            require(point.forceToStress(tangent,tangent,*segment).norm()<1.e-12,"Tangential force produces mobility stress");
            const auto resolved=point.forceToStress(glide,tangent,*segment);
            require((point.forceToStress(glide+100*tangent+20*normal,tangent,*segment)-resolved).norm()<1.e-10,
                    "Unresolved force changes mobility stress");
            const V recovered=(resolved*segment->burgers()).cross(tangent);
            require((recovered-glide).norm()<1.e-10,"Equivalent stress does not recover resolved force");
        }
        std::filesystem::remove_all(folder);
        std::cout << "PASS: quadrature glide direction and consistent resolved-force stress\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
