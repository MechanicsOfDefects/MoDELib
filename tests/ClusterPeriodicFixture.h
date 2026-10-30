#ifndef model_test_ClusterPeriodicFixture_H_
#define model_test_ClusterPeriodicFixture_H_
#include <ClusterDynamics.h>
#include <MicrostructureGenerator.h>
#include <DislocationNetwork.h>
#include <DislocationDynamicsBase.h>
#include <ConstrainedVelocitySolve.h>
#include <FixedDirichletSolver.h>
#include <filesystem>
#include <fstream>
#include <regex>
#include <chrono>
#include <iostream>
#include <iomanip>
using namespace model;
namespace fs = std::filesystem;
void require(bool ok, const std::string &message)
{
    if (!ok)
    {
        throw std::runtime_error(message);
    }
}
std::string read(const fs::path &path)
{
    std::ifstream stream(path);
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}
void parameter(std::string &text, const std::string &key, const std::string &value)
{
    const std::regex pattern("(^|\\n)\\s*" + key + "\\s*=[^;]*;");
    require(std::regex_search(text, pattern), "Missing parameter " + key);
    text = std::regex_replace(text, pattern, "$1" + key + "=" + value + ";");
}
fs::path fixture(const fs::path &root, int rank, int reactions, bool source, const std::string &name = "unitCube24.msh")
{
    const auto folder =
        fs::temp_directory_path() /
        ("cluster-periodic-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(folder / "inputFiles");
    fs::create_directories(folder / "config");
    fs::create_directories(folder / "F");
    fs::copy_file(root / "Library/Meshes" / name, folder / "inputFiles/mesh.msh");
    std::string material;
    if constexpr(ClusterDynamicsParameters<3>::mSize==1)
    {
        material=read(root / "Library/Materials/Zr_BMD19.txt");
    }
    else
    {
    material = read(root / "Library/Materials/Zr4_Fitted.txt");
    parameter(material, "computeReactions", reactions ? "1" : "0");
    parameter(material, "doseRate_dpaPerSec", source ? "1e-7" : "0");
    parameter(material, "otherSinks_SI", reactions == 1 ? "1e14 1e14 1e14 1e14" : "0 0 0 0");
    std::string map;
    const int species[4] = {-1, 1, 2, 3};
    for (int i = 0; i < 4; ++i)
    {
        for (int j = i; j < 4; ++j)
        {
            map += std::to_string(species[i]) + " " + std::to_string(species[j]) + " " +
                   ((reactions == 2 && i == 0 && j == 1) ? "1\n" : "0\n");
        }
    }
    parameter(material, "reactionPrefactorMap", map);
    if (reactions == 2)
    {
        parameter(material, "mobileSpeciesCascadeFractions", "1 1 0 0");
    }
    }
    std::ofstream(folder / "inputFiles/material.txt") << material;
    std::ofstream(folder / "inputFiles/ClusterDynamics.txt") << "useClusterDynamicsFEM=1;\n";
    const Eigen::Matrix3d F =
        (Eigen::Matrix3d() << 100, 50, 0, 0, 50 * std::sqrt(3.0), 0, 0, 0, 100 * std::sqrt(8.0 / 3.0)).finished();
    StaticID<Polyhedron<3, 3 - 1>>::force_count(0);
    std::set<int> seeds;
    {
        SimplicialMesh<3> mesh((root / "Library/Meshes" / name).string(), F, Eigen::Vector3d::Zero(), {-1});
        for (const auto &region : mesh.regions())
        {
            for (const auto &entry : region.second->faces())
            {
                if (entry.second->periodicFacePair.second)
                {
                    for (int d = 0; d < rank; ++d)
                    {
                        if (entry.second->periodicFacePair.first.normalized().cross(F.col(d).normalized()).norm() <
                            1.e-8)
                        {
                            seeds.insert(entry.first);
                        }
                    }
                }
            }
        }
    }
    StaticID<Polyhedron<3, 3 - 1>>::force_count(0);
    std::ofstream poly(folder / "inputFiles/polycrystal.txt");
    poly << std::setprecision(17) << "materialFile=material.txt;\nmeshFile=" << "mesh.msh"
         << ";\nabsoluteTemperature=600;\nF=";
    for (int i = 0; i < 3; ++i)
    {
        poly << F.row(i) << "\n";
    }
    poly << ";\nX0=0 0 0;\nperiodicFaceIDs=";
    for (int seed : seeds)
    {
        poly << seed << " ";
    }
    poly << ";\n";
    const int grains = name == "bicrystal_48.msh" ? 2 : (name.find("polycrystal") != std::string::npos ? 10 : 1);
    for (int g = 1; g <= grains; ++g)
    {
        poly << "C2G" << g << "=1 0 0 0 1 0 0 0 1;\n";
    }
    std::ofstream dc(folder / "inputFiles/DefectiveCrystal.txt");
    dc << "physics=ClusterDynamics;\nfmmTreeDepth=0;\nuseFEM=1;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;"
          "\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=";
    for (int d = 0; d < rank; ++d)
    {
        dc << "0 ";
    }
    dc << ";\nEwaldLengthFactor=0;\n";
    return folder;
}


#endif
