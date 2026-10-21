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
    auto material = read(root / "Library/Materials/Zr4_Fitted.txt");
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
    dc << "physics=ClusterDynamics;\nuseFEM=1;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;"
          "\noutputBinary=0;\nstartAtTimeStep=0;\nperiodicImageSize=";
    for (int d = 0; d < rank; ++d)
    {
        dc << "0 ";
    }
    dc << ";\nEwaldLengthFactor=0;\n";
    return folder;
}

void climb(const fs::path &root, int rank)
{
    const auto folder = fixture(root, rank, 0, false);
    auto dc = read(folder / "inputFiles/DefectiveCrystal.txt");
    parameter(dc, "physics", "ClusterDynamics DislocationDynamics");
    std::ofstream(folder / "inputFiles/DefectiveCrystal.txt") << dc;
    fs::copy_file(root / "Library/DislocationDynamics/DD.txt", folder / "inputFiles/DD.txt");
    DefectiveCrystalBase<3> base(folder.string());
    MicrostructureGenerator generator(base);
    const auto grain = base.poly.grain(1);
    const Eigen::Vector3d center = Eigen::Vector3d::Zero();
    const Eigen::Vector3d normal = Eigen::Vector3d::UnitZ();
    const auto key = GlidePlaneKey<3>(center, grain->reciprocalLatticeDirection(normal));
    const auto plane = base.dislocationDynamics().periodicGlidePlaneFactory.getFromKey(key);
    const Eigen::Vector3d burgers(0, 0, std::sqrt(8.0 / 3.0));
    const double halfWidth = rank > 0 ? 25.0 : 20.0;
    const double centerX = rank > 0 ? 40.0 : 0.0;
    std::vector<Eigen::Vector3d> polygon{{centerX - halfWidth, -10, 0},
                                         {centerX + halfWidth, -10, 0},
                                         {centerX + halfWidth, 10, 0},
                                         {centerX - halfWidth, 10, 0}};
    require(generator.insertJunctionLoop(polygon, plane, burgers, normal, center, 1, DislocationLoopIO<3>::SESSILELOOP),
            "Failed to insert climb loop");
    MicrostructureContainer<3> container(base);
    auto cd = std::make_shared<ClusterDynamics<3>>(container);
    container.push_back(cd);
    DislocationNetwork<3> network(container);
    network.initializeConfiguration(generator.config(), std::ofstream(), std::ofstream());
    const auto original = network.loops().begin()->second.lock();
    const auto duplicate = network.loops().create(original->burgers(), original->glidePlane);
    std::vector<std::shared_ptr<DislocationLoopNode<3>>> nodes;
    for (const auto *link : original->linkSequence())
    {
        const auto &node = link->source;
        nodes.push_back(network.loopNodes().create(duplicate, node->networkNode, node->get_P(),
                                                   node->periodicPlanePatch(), node->periodicPlaneEdge));
    }
    network.insertLoop(duplicate, nodes);
    network.updateGeometry();
    size_t quadrature = 0;
    for (const auto &entry : network.networkLinks())
    {
        const auto link = entry.second.lock();
        link->createQuadraturePoints(true);
        for (auto &qp : link->quadraturePoints())
        {
            qp.cCD.setZero();
            qp.cDD.setZero();
            qp.cDD(0) = 1.e-4;
            ++quadrature;
        }
    }
    require(quadrature > 0, "Climb regression has no active quadrature points");
    require(bool(network.climbSolver), "Climb solver not enabled");
    network.climbSolver->computeClimbScalarVelocities(true);
    const auto velocity = network.climbSolver->getNodeVelocities();
    require(velocity.allFinite() && velocity.norm() > 0, "Climb drive produced no finite motion");
    size_t relations = 0;
    forEachPeriodicVelocityRelation(
        network,
        [&](size_t n, size_t a, size_t b, double wa, double wb)
        {
            ++relations;
            require((velocity.segment<3>(3 * n) - wa * velocity.segment<3>(3 * a) - wb * velocity.segment<3>(3 * b))
                            .norm() < 1.e-9 * velocity.norm(),
                    "Actual climb solver violates periodic vector interpolation");
        });
    require(rank == 0 || relations > 0, "Periodic climb test has no crossings");
    size_t tangentialMarkers = 0;
    size_t localIndex = 0;
    for (const auto &entry : network.networkNodes())
    {
        const auto node = entry.second.lock();
        const auto direction = node->climbDirection();
        const auto value = velocity.segment<3>(3 * localIndex);
        const auto &species = network.climbSolver->scalarVelocities()[localIndex];
        const double normalSpeed = -(species * cd->cdp.msVector / cd->cdp.msVector.abs()).sum();
        require(std::abs(value.dot(direction) - normalSpeed) < 1.e-9 * velocity.norm(),
                "Climb source strength differs from normal motion");
        if ((value - direction * value.dot(direction)).norm() > 1.e-8 * velocity.norm())
        {
            ++tangentialMarkers;
        }
        ++localIndex;
    }
    require(rank == 0 || tangentialMarkers > 0, "Regression does not exercise tangential subdivision motion");
    network.climbSolver->computeClimbScalarVelocities(false);
    require((network.climbSolver->getNodeVelocities() - velocity).norm() < 1.e-10 * velocity.norm(),
            "Cached climb solve changed velocities");
    size_t index = 0;
    for (const auto &entry : network.networkNodes())
    {
        entry.second.lock()->set_V(velocity.segment<3>(3 * index++), true);
    }
    network.moveNodes(1.e-4 / velocity.norm());
    for (const auto &entry : network.networkLinks())
    {
        require(entry.second.lock()->loopLinks().size() == 2, "Climb motion broke shared junction ownership");
    }
    fs::remove_all(folder);
}
int main(int argc, char **argv)
{
    try
    {
        require(argc == 2, "Expected repository root");
        for (int rank : {0, 1, 3})
        {
            climb(argv[1], rank);
        }
        std::cout << "PASS: coupled CD/climb solver, periodic crossings, motion and junction sharing\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
