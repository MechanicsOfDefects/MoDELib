#include <EllipsoidalInclusionIndividualSpecification.h>
#include <EllipsoidalInclusion.h>
#include <filesystem>
#include <chrono>
#include <fstream>
#include <stdexcept>

void check(bool ok, const char *message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
int main(int argc, char **argv)
{
    try
    {
        using Spec = model::EllipsoidalInclusionIndividualSpecification;
        check(argc == 2, "Expected example filename");
        const Spec defaults;
        check(defaults.type == "EllipsoidalInclusion" && defaults.style == "Individual" && !defaults.parser &&
                  defaults.center.isZero() && defaults.crystalAxes_SI.isZero() && defaults.crystalEigenDistortion.isZero() &&
                  defaults.velocityReductionFactor == 1.0 && defaults.phaseID == 0 && !defaults.allowOverlap &&
                  !defaults.allowOutside,
              "Incorrect specification defaults");
        const Spec parsed(argv[1]);
        Eigen::Matrix3d expectedAxes;
        expectedAxes << 0, -2e-8, 0, 3e-8, 0, 0, 0, 0, 1e-8;
        Eigen::Matrix3d expectedStrain;
        expectedStrain << 1e-3, 2e-4, 0, 2e-4, -2e-3, 0, 0, 0, 3e-3;
        check(parsed.center.isApprox(Eigen::Vector3d(200, 0, 0)) && parsed.crystalAxes_SI.isApprox(expectedAxes) &&
                  parsed.crystalEigenDistortion(1) == 2e-4 && parsed.crystalEigenDistortion(8) == 3e-3 &&
                  parsed.velocityReductionFactor == 1.0 && parsed.phaseID == 0 && !parsed.allowOverlap &&
                  !parsed.allowOutside,
              "Example parsing failed");
        const Eigen::Matrix3d strain = Eigen::Map<const Eigen::Matrix3d>(parsed.crystalEigenDistortion.data()).transpose();
        check(strain.isApprox(expectedStrain), "Eigen-distortion ordering failed");
        const double b_SI = 2.5e-10;
        const model::EllipsoidalInclusion<3> inclusion(parsed.center, parsed.crystalAxes_SI / b_SI, strain, .31, 43.0,
                                                       parsed.velocityReductionFactor, parsed.phaseID, nullptr);
        check(inclusion.contains(parsed.center) && inclusion.stress(parsed.center).allFinite(),
              "Specification cannot construct an inclusion");

        const auto file = std::filesystem::temp_directory_path() /
                          ("model-ellipsoid-spec-" +
                           std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".txt");
        {
            std::ofstream out(file);
            out << "center=1 2 3;\ncrystalAxes_SI=1 2;\n";
        }
        bool rejected = false;
        try
        {
            const Spec invalid(file.string());
        }
        catch(const std::runtime_error &)
        {
            rejected = true;
        }
        check(rejected, "Wrong axes dimensions accepted");
        {
            std::ofstream out(file);
            out << "center=1 2 3;\ncrystalAxes_SI=1e-8 0 0 0 2e-8 0 0 0 3e-8;\n"
                << "crystalEigenDistortion=0 0 0 0 0 0 0 0 0;\n"
                << "velocityReductionFactor=.75;\nphaseID=4;\nallowOverlap=1;\nallowOutside=1;\n";
        }
        const Spec custom(file.string());
        check(custom.center.isApprox(Eigen::Vector3d(1, 2, 3)) && custom.velocityReductionFactor == .75 &&
                  custom.phaseID == 4 && custom.allowOverlap && custom.allowOutside,
              "Custom metadata parsing failed");
        std::filesystem::remove(file);
        std::cout
            << "PASS specification defaults, file parsing, field order, unit conversion, and malformed dimensions\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
