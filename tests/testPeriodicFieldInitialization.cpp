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
            ("periodic-fields-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
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
        for(int rank = 0; rank <= 3; ++rank)
        {
            for(int imageCount : {0, 1})
            {
                {
                    std::ofstream f(input / "DefectiveCrystal.txt");
                    f << "physics=ElasticDeformation;\nuseFEM=0;\nNsteps=0;\nmaxResolveSteps=0;\ndtMax=1;\n"
                      << "outputFrequency=1;\noutputBinary=0;\nstartAtTimeStep=0;\n";
                    if(rank > 0)
                    {
                        f << "periodicImageSize=";
                        for(int i = 0; i < rank; ++i)
                        {
                            f << imageCount << ' ';
                        }
                        f << ";\nEwaldLengthFactor=0.5;\n";
                    }
                    std::ofstream p(input / "polycrystal.txt");
                    const double length = 200.0 / std::sqrt(3.0);
                    p << std::setprecision(17) << "materialFile=W.txt;\nmeshFile=" << "mesh.msh"
                      << ";\nabsoluteTemperature=300;\nC2G1=1 0 0 0 1 0 0 0 1;\nF="
                      << length << " 0 0 0 " << length << " 0 0 0 " << length << ";\nX0=0 0 0;\nperiodicFaceIDs=";
                    for(int i = 0; i < rank; ++i)
                    {
                        p << seeds[i] << ' ';
                    }
                    p << ";\n";
                }
                model::DislocationDynamicsBase<3> base(folder.string());
                require(base.periodicLatticeBasis.cols() == rank, "Incorrect active rank");
                require(base.isPeriodicDomain == (rank == 3), "Legacy flag changed");
                require(base.periodicImageSize.size() == size_t(rank), "Image counts ignored for active directions");
                const size_t expected = imageCount == 0 ? 1 : size_t(std::pow(3, rank));
                require(base.periodicShifts.size() == expected, "Incorrect image count");
                std::set<std::vector<int>> coefficients;
                for(const auto& shift : base.periodicShifts)
                {
                    const Eigen::VectorXd coordinates = base.periodicLatticeReciprocalBasis.transpose() * shift;
                    std::vector<int> key;
                    for(int i = 0; i < rank; ++i)
                    {
                        key.push_back(std::lround(coordinates[i]));
                        require(std::abs(coordinates[i] - key.back()) < 1.e-8 && std::abs(key.back()) <= imageCount,
                                "Invalid periodic shift");
                    }
                    coefficients.insert(key);
                }
                require(coefficients.size() == expected, "Duplicate periodic shift");
                require(rank == 0 ? base.EwaldLength == 0.0 : base.EwaldLength > 0.0,
                        "Ewald initialization depends on domain flag or image count");
            }
        }
        std::filesystem::remove_all(folder);
        std::cout << "PASS: image and Ewald initialization for zero, one, two and three periodic directions\n";
        return 0;
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
