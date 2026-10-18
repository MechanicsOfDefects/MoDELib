#include <MicrostructureGenerator.h>
#include <DislocationNetwork.h>
#include <MicrostructureContainer.h>
#include <DislocationDynamicsBase.h>
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
        model::DislocationDynamicsBase<3> base(folder.string());
        model::MicrostructureGenerator generator(base);
        model::PrismaticLoopIndividualSpecification spec;
        spec.slipSystemIDs={0};
        spec.loopRadii={1.e-7};
        spec.glideSteps={50};
        spec.loopCenters.resize(1,3);
        spec.loopCenters << -0.5*length,0,0;
        generator.addPrismaticLoopIndividual(spec);
        require(generator.config().loops().size()==7,"Expected complete prismatic construction");
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        network.initializeConfiguration(generator.config(),std::ofstream(),std::ofstream());
        network.alphaLineTension=1.0;
        require(!network.stochasticForceGenerator,"Regression requires deterministic mobility");
        using V=Eigen::Vector3d;
        size_t checked=0;
        double maxRawTangential=0;
        for(const auto& entry : network.networkLinks())
        {
            const auto segment=entry.second.lock();
            segment->createQuadraturePoints(false);
            segment->updateQuadraturePoints(false);
            if(!segment->slipSystem())
            {
                continue;
            }
            const V normal=segment->slipSystem()->unitNormal;
            for(auto& point : segment->quadraturePoints())
            {
                ++checked;
                const V glide=normal.cross(point.rl).normalized();
                const V raw=point.lineTensionForce;
                maxRawTangential=std::max(maxRawTangential,std::abs(raw.dot(point.rl)));
                require(std::abs(point.pkForce.dot(glide))<1.e-9,"Initial prism has unexpected PK glide force");
                require(std::abs(raw.dot(glide))<1.e-9,"Initial prism has unexpected LT glide force");
                require(point.velocity.norm()<1.e-12,"Unloaded bicrystal prism has spurious glide velocity");
                // Confirm raw line tension still equals energy times spline curvature.
                V expected=V::Zero();
                for(const auto* link : segment->loopLinks())
                {
                    const double u=link->source->networkNode==segment->source ? point.abscissa : 1-point.abscissa;
                    const V curvature=link->spline().get_rll(u);
                    if(!std::isnan(curvature.squaredNorm()))
                    {
                        expected+=point.coreEnergyPerLength*curvature;
                    }
                }
                require((raw-expected).norm()<1.e-12,"Raw line tension was altered");
                const auto initialStress=point.internalStress;
                // Apply an actual shear stress producing a specified PK glide force.
                const auto drive=point.forceToStress(1.e-3*glide,point.rl,*segment);
                point.internalStress=initialStress+drive;
                point.updateForcesAndVelocities(*segment,false);
                const V positive=point.velocity;
                require(positive.norm()>1.e-8 && positive.dot(glide)>0,"Positive glide force gives no forward response");
                require(std::abs(positive.dot(point.rl))<1.e-12 && std::abs(positive.dot(normal))<1.e-12,
                        "Driven velocity is not along the permitted glide direction");
                require((point.lineTensionForce-raw).norm()<1.e-12,"Velocity update modified raw line tension");
                // Unresolved line tension must change neither speed nor direction.
                point.lineTensionForce.setZero();
                point.updateForcesAndVelocities(*segment,false);
                require((point.velocity-positive).norm()<1.e-9,"Climb/tangential LT changes driven mobility");
                point.lineTensionForce=raw;
                point.internalStress=initialStress-drive;
                point.updateForcesAndVelocities(*segment,false);
                require(point.velocity.norm()>1.e-8 && point.velocity.dot(glide)<0,"Negative glide force gives no reverse response");
                require(std::abs(point.velocity.dot(point.rl))<1.e-12,"Reverse velocity is tangential");
                point.internalStress=initialStress;
                point.updateForcesAndVelocities(*segment,false);
            }
        }
        require(checked>0,"No quadrature velocities tested");
        require(maxRawTangential>1.e-6,"Test did not reproduce a finite tangential LT component");
        std::filesystem::remove_all(folder);
        std::cout << "PASS: " << checked << " bicrystal prism quadrature points, unloaded and driven, raw LT unchanged\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
