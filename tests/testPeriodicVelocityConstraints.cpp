#include <DislocationNetwork.h>
#include <PeriodicVelocityConstraints.h>
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
            ("bicrystal-shared-constraints-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        model::MicrostructureContainer<3> container(base);
        model::DislocationNetwork<3> network(container);
        model::DDconfigIO<3> saved(folder.string());
        std::ifstream stream(root/"tests/data/bicrystalSharedConstraints1915.txt");
        require(bool(stream),"Missing bicrystal junction regression fixture");
        saved.readTxtStream(stream);
        network.initializeConfiguration(saved,std::ofstream(),std::ofstream());
        using Triplets=std::vector<Eigen::Triplet<double>>;
        Triplets triplets;
        const size_t freeCount=model::periodicVelocityConstraints<3>(network,triplets);
        Eigen::SparseMatrix<double> sparse(3*network.networkNodes().size(),3*freeCount);
        sparse.setFromTriplets(triplets.begin(),triplets.end());
        const Eigen::MatrixXd Z(sparse);
        require((Z*Eigen::VectorXd::Ones(Z.cols())-Eigen::VectorXd::Ones(Z.rows())).norm()<1.e-8,
                "Periodic constraints lost uniform translation or doubled a velocity");
        size_t shared=0;
        for(const auto& entry : network.networkNodes())
        {
            const auto node=entry.second.lock();
            if(node->sID>=975875 && node->sID<=975878)
            {
                ++shared;
                require(node->loopNodes().size()==2,"Missing shared crossing memberships");
                for(const auto* ln : node->loopNodes())
                {
                    require(ln->isPeriodicSubdivision(),"Fixture no longer exercises periodic interpolation");
                    const auto a=ln->periodicPrev();
                    const auto b=ln->periodicNext();
                    const double da=(ln->get_P()-a->get_P()).norm();
                    const double db=(ln->get_P()-b->get_P()).norm();
                    const auto residual=(Z.middleRows(3*node->networkID(),3)
                        -db/(da+db)*Z.middleRows(3*a->networkNode->networkID(),3)
                        -da/(da+db)*Z.middleRows(3*b->networkNode->networkID(),3)).eval();
                    require(residual.norm()<1.e-8,"Shared crossing does not satisfy both loop interpolants");
                }
            }
        }
        require(shared==4,"Missing frame 1915 shared crossings");
        // Independent endpoints 1 and 2 must become compatible through a
        // shared crossing. Repeated/reversed relations must not add rank.
        for(bool reverse : {false,true})
        {
            model::PeriodicVelocityConstraints constraints(5);
            for(size_t endpoint : {reverse ? 2UL : 1UL,reverse ? 1UL : 2UL})
            {
                constraints.interpolate(3,0,endpoint,0.4,0.6);
                constraints.interpolate(4,endpoint,0,0.2,0.8);
                constraints.interpolate(3,endpoint,0,0.6,0.4);
            }
            Triplets t;
            const auto count=constraints.append<1>(t);
            require(count==2,"Duplicate periodic relations changed the number of free velocities");
            Eigen::SparseMatrix<double> transform(5,count);
            transform.setFromTriplets(t.begin(),t.end());
            const Eigen::MatrixXd M(transform);
            require((M.row(1)-M.row(2)).norm()<1.e-12,"Competing endpoint velocities are incompatible");
            require((M.row(3)-0.4*M.row(0)-0.6*M.row(1)).norm()<1.e-12,"Incorrect crossing interpolation");
            require((M.row(4)-0.8*M.row(0)-0.2*M.row(2)).norm()<1.e-12,"Incorrect second crossing interpolation");
            require((M*Eigen::Vector2d::Ones()-Eigen::VectorXd::Ones(5)).norm()<1.e-12,"Constant velocity changed");
        }
        std::filesystem::remove_all(folder);
        std::cout<<"PASS: frame 1915 shared periodic constraints, compatibility and uniform translation\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
