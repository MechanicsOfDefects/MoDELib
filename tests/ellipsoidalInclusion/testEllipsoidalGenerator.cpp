#include <MicrostructureGenerator.h>
#include <InclusionMicrostructure.h>
#include <EllipsoidalInclusion.h>
#include <filesystem>
#include <fstream>
#include <chrono>
#include <random>
#include <stdexcept>
using V = Eigen::Vector3d;
using M = Eigen::Matrix3d;
void check(bool ok, const char *msg)
{
    if(!ok)
    {
        throw std::runtime_error(msg);
    }
}
int main(int argc, char **argv)
{
    try
    {
        check(argc == 2, "Expected repository root");
        const std::filesystem::path root(argv[1]);
        const auto dir = std::filesystem::temp_directory_path() /
                         ("model-ellipsoid-generator-" +
                          std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        for(const auto sub : {"inputFiles", "evl", "aux", "F"})
        {
            std::filesystem::create_directories(dir / sub);
        }
        const auto input = dir / "inputFiles";
        std::filesystem::copy_file(root / "Library/Materials/W.txt", input / "W.txt");
        std::filesystem::copy_file(root / "Library/Meshes/unitCube24.msh", input / "mesh.msh");
        {
            std::ofstream f(input / "DefectiveCrystal.txt");
            f << "physics=none;\nuseFEM=0;\nNsteps=1;\nmaxResolveSteps=0;\ndtMax=1;\noutputFrequency=1;\noutputBinary="
                 "0;\nstartAtTimeStep=0;\n";
        }
        {
            std::ofstream f(input / "polycrystal.txt");
            f << "materialFile=W.txt;\nmeshFile=mesh.msh;\nabsoluteTemperature=300;\n"
              << "C2G1=0.6 -0.8 0 0.8 0.6 0 0 0 1;\nF=1000 0 0 0 1000 0 0 0 1000;\nX0=0 0 0;\nperiodicFaceIDs=;\n";
        }
        model::DislocationDynamicsBase<3> base(dir.string());
        model::MicrostructureGenerator mg(base);
        model::EllipsoidalInclusionIndividualSpecification spec;
        spec.center = base.mesh.xCenter();
        const M axes = V(30, 20, 10).asDiagonal();
        M C2G;
        C2G << .6, -.8, 0, .8, .6, 0, 0, 0, 1;
        const M globalAxes = C2G * axes;
        spec.crystalAxes_SI = axes * base.poly.b_SI;
        spec.crystalEigenDistortion << .01, .002, 0, .002, .01, 0, 0, 0, .01;
        const M crystalStrain = Eigen::Map<const M>(spec.crystalEigenDistortion.data()).transpose();
        const M globalStrain = C2G * crystalStrain * C2G.transpose();
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 1, "Valid inclusion rejected");
        check(mg.config().ellipsoidalInclusions()[0].eT.isApprox(globalStrain), "Ellipsoid tensor rotation failed");
        check(mg.config().ellipsoidalInclusions()[0].axes.isApprox(globalAxes),
              "Wrong crystal-to-global rotation or unit conversion");
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 1, "Overlap accepted");
        spec.allowOverlap = true;
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 2, "allowOverlap ignored");
        spec.center = base.mesh.xMax() + V::Ones() * 1000;
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 2, "Outside inclusion accepted");
        spec.allowOutside = true;
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 2, "Outside center has no crystal frame");
        spec.center = V(490, 0, 0);
        spec.allowOutside = false;
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 2, "Boundary crossing accepted");
        spec.allowOutside = true;
        mg.addEllipsoidalInclusionIndividual(spec);
        check(mg.config().ellipsoidalInclusions().size() == 3, "allowOutside ignored for inside center");
        spec.crystalAxes_SI.setZero();
        bool rejected = false;
        try
        {
            mg.addEllipsoidalInclusionIndividual(spec);
        }
        catch(const std::invalid_argument &)
        {
            rejected = true;
        }
        check(rejected && mg.config().ellipsoidalInclusions().size() == 3, "Invalid input mutated configuration");
        // Density generation: independent scales, crystal-frame planes, and target count.
        model::EllipsoidalInclusionDensitySpecification cloud;
        model::MicrostructureGenerator cloudGenerator(base);
        cloudGenerator.addEllipsoidalInclusionDensity(cloud);
        check(cloudGenerator.config().ellipsoidalInclusions().empty(), "Zero density generated inclusions");
        const double volumeSI = base.mesh.volume() * std::pow(base.poly.b_SI, 3);
        cloud.targetDensity = 11.5 / volumeSI;
        cloud.crystalAxes = M::Identity();
        cloud.crystalEigenDistortion = spec.crystalEigenDistortion;
        cloud.axesLognormalDistributionMean_SI = V(30, 24, 8).transpose() * base.poly.b_SI;
        cloud.axesLognormalDistributionStd_SI = V(3, 4.8, 2.4).transpose() * base.poly.b_SI;
        cloud.patternVector_SI << 100 * base.poly.b_SI, 0, 0;
        cloud.allowOutside = true;
        cloud.allowOverlap = true;
        cloudGenerator.addEllipsoidalInclusionDensity(cloud);
        check(cloudGenerator.config().ellipsoidalInclusions().size() == 12, "Density count incorrect");
        const V planeNormal = C2G.col(0);
        std::mt19937 expectedRandom;
        std::normal_distribution<double> expectedNormal(0.0, 1.0);
        bool independentScales = false;
        V firstLengths = cloudGenerator.config().ellipsoidalInclusions().front().axes.colwise().norm();
        bool variedSizes = false;
        for(const auto &inclusion : cloudGenerator.config().ellipsoidalInclusions())
        {
            const double planeIndex = inclusion.C.dot(planeNormal) / 100;
            check(std::abs(planeIndex - std::round(planeIndex)) < 1.e-10, "Center missed crystal-frame plane");
            const M localAxes = C2G.transpose() * inclusion.axes;
            const V lengths = localAxes.colwise().norm();
            check((localAxes * lengths.cwiseInverse().asDiagonal()).isApprox(M::Identity(), 1.e-12),
                  "Density generation rotated crystal semiaxes incorrectly");
            check(inclusion.eT.isApprox(globalStrain), "Density tensor rotation failed");
            for(int k = 0; k < 3; ++k)
            {
                const double mean = cloud.axesLognormalDistributionMean_SI(k);
                const double cv = cloud.axesLognormalDistributionStd_SI(k) / mean;
                const double variance = std::log1p(cv * cv);
                const double expectedLength =
                    std::exp(std::log(mean) - variance / 2 + std::sqrt(variance) * expectedNormal(expectedRandom));
                check(std::abs(lengths(k) * base.poly.b_SI / expectedLength - 1) < 1.e-12,
                      "Arithmetic SI mean/std were not converted to lognormal parameters");
            }
            const V factors = lengths.cwiseQuotient(V(30, 20, 10));
            independentScales = independentScales || std::abs(factors(0) - factors(1)) > 1.e-6;
            variedSizes = variedSizes || !lengths.isApprox(firstLengths);
        }
        check(independentScales && variedSizes, "Semiaxes were not perturbed independently");
        cloud.targetDensity = .5 / volumeSI;
        cloud.axesLognormalDistributionStd_SI.setZero();
        cloud.patternVector_SI.setZero();
        cloudGenerator.addEllipsoidalInclusionDensity(cloud);
        check(cloudGenerator.config().ellipsoidalInclusions().size() == 13, "Density must add new inclusions");
        check(cloudGenerator.config().ellipsoidalInclusions().back().axes.isApprox(
                  C2G * (cloud.axesLognormalDistributionMean_SI / base.poly.b_SI).asDiagonal()),
              "Zero-spread sizes incorrect");
        cloud.axesLognormalDistributionStd_SI(0) = -1;
        bool invalidDensity = false;
        try
        {
            cloudGenerator.addEllipsoidalInclusionDensity(cloud);
        }
        catch(const std::invalid_argument &)
        {
            invalidDensity = true;
        }
        check(invalidDensity && cloudGenerator.config().ellipsoidalInclusions().size() == 13,
              "Invalid density distribution changed configuration");
        cloud.axesLognormalDistributionStd_SI.setZero();
        cloud.axesLognormalDistributionMean_SI.setConstant(2000 * base.poly.b_SI);
        cloud.allowOutside = false;
        cloud.maxAttempts = 2;
        bool exhausted = false;
        try
        {
            cloudGenerator.addEllipsoidalInclusionDensity(cloud);
        }
        catch(const std::runtime_error &)
        {
            exhausted = true;
        }
        check(exhausted && cloudGenerator.config().ellipsoidalInclusions().size() == 13,
              "Impossible density placement did not stop");
        // File parsing and density dispatch, including nonzero default spread/pattern.
        {
            std::ofstream f(input / "ellipsoidDensity.txt");
            f << std::setprecision(17) << "type=EllipsoidalInclusions;\nstyle=density;\ntargetDensity=" << .5 / volumeSI
              << ";\n"
              << "crystalAxes=" << M::Identity() << ";\n"
              << "axesLognormalDistributionMean_SI=" << V(30, 20, 10).transpose() * base.poly.b_SI
              << ";\naxesLognormalDistributionStd_SI=0 0 0;\n"
              << "crystalEigenDistortion=.01 .002 0 .002 .01 0 0 0 .01;\n"
              << "patternVector_SI=0 0 0;\nallowOverlap=1;\nallowOutside=1;\n"
              << "velocityReductionFactor=.75;\nphaseID=0;\nmaxAttempts=100;\n";
            std::ofstream list(input / "initialMicrostructure.txt");
            list << "microstructureFile=ellipsoidDensity.txt;\n";
        }
        cloudGenerator.readMicrostructureFile();
        check(cloudGenerator.config().ellipsoidalInclusions().size() == 1 &&
                  cloudGenerator.config().ellipsoidalInclusions().front().axes.isApprox(globalAxes) &&
                  cloudGenerator.config().ellipsoidalInclusions().front().eT.isApprox(globalStrain),
              "Density file dispatch failed");
        const model::EllipsoidalInclusionDensitySpecification templateDensity(
            (root / "Library/Microstructures/ellipsoidalInclusionsDensity.txt").string());
        check(templateDensity.axesLognormalDistributionStd_SI(0) == 6e-9 && templateDensity.maxAttempts == 10000,
              "Density library template parsing failed");
        // Exercise both spherical entry points and polyhedral generation in a rotated grain.
        model::MicrostructureGenerator other(base);
        model::SphericalInclusionIndividualSpecification spherical;
        spherical.centers = V::Zero().transpose();
        spherical.radii_SI = {5 * base.poly.b_SI};
        spherical.crystalEigenDistortions = spec.crystalEigenDistortion;
        spherical.velocityReductionFactors = {1.0};
        spherical.phaseIDs = {0};
        spherical.allowOverlap = true;
        spherical.allowOutside = false;
        other.addSphericalInclusionIndividual(spherical);
        check(other.config().sphericalInclusions().size() == 1 &&
                  other.config().sphericalInclusions().back().eT.isApprox(globalStrain),
              "Individual sphere tensor rotation failed");
        model::SphericalInclusionDensitySpecification density;
        density.targetDensity = 0.5 / (base.mesh.volume() * std::pow(base.poly.b_SI, 3));
        density.diameterLognormalDistribution_M = 2 * base.poly.b_SI;
        density.diameterLognormalDistribution_A = base.poly.b_SI;
        density.diameterLognormalDistribution_S = 0;
        density.crystalEigenDistortion = spec.crystalEigenDistortion;
        density.allowOverlap = true;
        density.allowOutside = true;
        density.phaseID = 0;
        other.addSphericalInclusionDensity(density);
        check(other.config().sphericalInclusions().size() == 2 &&
                  other.config().sphericalInclusions().back().eT.isApprox(globalStrain),
              "Density sphere tensor rotation failed");
        model::PolyhedronInclusionIndividualSpecification polyhedral;
        polyhedral.mshFile = (input / "mesh.msh").string();
        polyhedral.X0 = V::Zero().transpose();
        polyhedral.F = M::Identity() * 10;
        polyhedral.crystalEigenDistortion = spec.crystalEigenDistortion;
        polyhedral.velocityReductionFactor = 1;
        polyhedral.phaseID = 0;
        other.addPolyhedronInclusionIndividual(polyhedral);
        check(other.config().polyhedronInclusions().size() == 1 &&
                  other.config().polyhedronInclusions().back().eT.isApprox(globalStrain),
              "Polyhedron tensor rotation failed");
        const size_t sphere = mg.insertInclusion(V::Zero(), 1., M::Zero(), 1., 0);
        check(sphere == 3, "Sphere ID collided with ellipsoid");
        mg.config().ellipsoidalInclusions()[0].inclusionID = 100;
        check(mg.insertInclusion(V::Zero(), axes, M::Zero(), 1., 0) == 101, "Sparse IDs collided");
        std::map<size_t, V> nodes{{0, V::Zero()}, {1, V::UnitX()}, {2, V::UnitY()}};
        std::map<size_t, std::vector<size_t>> faces{{0, {0, 1, 2}}};
        check(mg.insertInclusion(nodes, faces, M::Zero(), 1., 0) == 102, "Polyhedron ID collided");
        // File dispatcher, followed by a round trip through generated records.
        {
            std::ofstream f(input / "ellipsoid.txt");
            f << std::setprecision(17)
              << "type=EllipsoidalInclusions;\nstyle=individual;\ncenter=" << base.mesh.xCenter().transpose() << ";\n"
              << "crystalAxes_SI=" << axes * base.poly.b_SI << ";\n"
              << "crystalEigenDistortion=.01 .002 0 .002 .01 0 0 0 "
                 ".01;\nvelocityReductionFactor=.75;\nphaseID=0;\nallowOverlap=0;\nallowOutside=0;\n";
            std::ofstream list(input / "initialMicrostructure.txt");
            list << "microstructureFile=ellipsoid.txt;\n";
        }
        mg.readMicrostructureFile();
        check(mg.config().ellipsoidalInclusions().size() == 1, "File dispatch failed");
        check(mg.config().ellipsoidalInclusions()[0].axes.isApprox(globalAxes), "File input lost grain rotation");
        check(mg.config().ellipsoidalInclusions()[0].mobilityReduction == .75, "Metadata lost");
        mg.config().writeTxt(0);
        model::DDconfigIO<3> restored((dir / "evl").string());
        restored.readTxt(0);
        check(restored.ellipsoidalInclusions().size() == 1, "Generated inclusion not serialized");
        check(restored.ellipsoidalInclusions()[0].eT.isApprox(globalStrain), "Serialized tensor rotation failed");
        // Load generated records into the simulation and compare its stress sum.
        model::MicrostructureContainer<3> simulation(base);
        auto inclusions = std::make_shared<model::InclusionMicrostructure<3>>(simulation);
        simulation.microstructures().push_back(inclusions);
        std::ofstream unused;
        // A sphere with a larger ID is read first, testing interleaved shape IDs.
        restored.sphericalInclusions().emplace_back(7, base.mesh.xCenter() + V(100, 0, 0), 5.0, M::Identity() * .002,
                                                    1.0, 0);
        simulation.initializeConfiguration(restored, unused, unused);
        check(inclusions->eshelbyInclusions().size() == 2, "Simulation failed to load mixed inclusions");
        const auto &record = restored.ellipsoidalInclusions().front();
        const auto loaded = std::dynamic_pointer_cast<model::EllipsoidalInclusion<3>>(
            inclusions->eshelbyInclusions().at(record.inclusionID));
        check(loaded && loaded->sID == record.inclusionID && loaded->phaseID == record.phaseID &&
                  loaded->mobilityReduction == record.mobilityReduction && loaded->secondPhase,
              "Loaded ellipsoid metadata incorrect");
        const model::EllipsoidalInclusion<3> reference(record.C, record.axes, record.eT, base.poly.nu, base.poly.mu,
                                                       record.mobilityReduction, record.phaseID, loaded->secondPhase);
        const auto &sphereRecord = restored.sphericalInclusions().front();
        const model::SphericalInclusion<3> sphereReference(sphereRecord.C, sphereRecord.a, sphereRecord.eT,
                                                           base.poly.nu, base.poly.mu, 1.0, 0, loaded->secondPhase);
        for(const V &point : {record.C, V(record.C + V(70, 45, 20))})
        {
            M expected = M::Zero();
            for(const auto &shift : base.periodicShifts)
            {
                expected += reference.stress(point + shift) + sphereReference.stress(point + shift);
            }
            check(expected.norm() > 0 && simulation.stress(point, nullptr, nullptr, nullptr).isApprox(expected, 1.e-12),
                  "Ellipsoid stress missing from simulation sum");
        }
        model::DDconfigIO<3> saved((dir / "evl").string());
        model::DDauxIO<3> aux((dir / "aux").string());
        simulation.output(saved, aux, unused, unused);
        check(saved.ellipsoidalInclusions().size() == 1 && saved.sphericalInclusions().size() == 1,
              "Simulation output lost an inclusion");
        check(saved.ellipsoidalInclusions()[0].axes.isApprox(record.axes) &&
                  saved.ellipsoidalInclusions()[0].eT.isApprox(record.eT) &&
                  saved.ellipsoidalInclusions()[0].inclusionID == record.inclusionID,
              "Simulation output changed ellipsoid data");
        saved.writeBin(1);
        model::DDconfigIO<3> restart((dir / "evl").string());
        restart.readBin(1);
        const V probe = record.C + V(70, 45, 20);
        const M before = simulation.stress(probe, nullptr, nullptr, nullptr);
        simulation.initializeConfiguration(restart, unused, unused);
        check(simulation.stress(probe, nullptr, nullptr, nullptr).isApprox(before, 1.e-12),
              "Restart changed the stress field");
        for(const int phaseID : {-1, 999})
        {
            restart.ellipsoidalInclusions()[0].phaseID = phaseID;
            bool invalidPhase = false;
            try
            {
                simulation.initializeConfiguration(restart, unused, unused);
            }
            catch(const std::runtime_error &)
            {
                invalidPhase = true;
            }
            check(invalidPhase, "Invalid ellipsoidal phase accepted");
        }
        restart.ellipsoidalInclusions()[0].phaseID = 0;
        restart.ellipsoidalInclusions()[0].inclusionID = 7;
        bool duplicate = false;
        try
        {
            simulation.initializeConfiguration(restart, unused, unused);
        }
        catch(const std::runtime_error &)
        {
            duplicate = true;
        }
        check(duplicate, "Duplicate inclusion IDs accepted");
        restart.clear();
        simulation.initializeConfiguration(restart, unused, unused);
        check(inclusions->eshelbyInclusions().empty() && simulation.stress(probe, nullptr, nullptr, nullptr).isZero(),
              "Empty reload retained stress");
        std::filesystem::remove_all(dir);
        std::cout
            << "PASS generator insertion, IDs, flags, validation, file dispatch, simulation stress and restart IO\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
