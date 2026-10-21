#include <DefectiveCrystal/DefectiveCrystalConfigIO.h>
#include <array>
#include <filesystem>
#include <sstream>
#include <stdexcept>

using Config = model::DefectiveCrystalConfigIO<3>;
void check(bool ok, const char *message)
{
    if(!ok)
    {
        throw std::runtime_error(message);
    }
}
std::string text(Config &config)
{
    std::ostringstream out;
    config.writeTxtStream(out);
    return out.str();
}
int main()
{
    try
    {
        Config defaultConfig;
        check(defaultConfig.getTxtFilename(7) == "config/evl_7.txt", "Default configuration folder is wrong");
        const auto directory =
            std::filesystem::temp_directory_path() /
            ("model-ddconfig-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directory(directory);
        Config config(directory.string()), restored(directory.string());
        Eigen::Matrix3d axes;
        axes << 0, -2, 0, 3, 0, 0, 0, 0, 1;
        Eigen::Matrix3d strain;
        strain << .01, .003, -.004, .003, -.002, .005, -.004, .005, .007;
        const Eigen::Vector3d center(1, -2, 3);
        config.sphericalInclusions().emplace_back(10, center, 2.0, strain, .5, 1);
        config.ellipsoidalInclusions().emplace_back(20, center, axes, strain, .625, 2);
        config.ellipsoidalInclusions().emplace_back(21, -center, 2 * axes, -strain, .75, 3);
        for(size_t id = 30; id < 33; ++id)
        {
            config.polyhedronInclusions().emplace_back(id, strain, .875, 4);
        }
        config.polyhedronInclusionNodes().emplace_back(100, center);
        config.polyhedronInclusionEdges().emplace_back(30, 0, 100, 100);
        const std::array<size_t, 11> counts{{0, 0, 0, 0, 1, 2, 3, 1, 1, 0, 0}};
        const auto serialized = text(config);
        std::istringstream header(serialized);
        std::string line;
        for(auto count : counts)
        {
            std::getline(header, line);
            check(std::stoull(line) == count, "Text count order is wrong");
        }
        for(auto id : {10, 20, 21, 30, 31, 32, 100, 30})
        {
            std::getline(header, line);
            check(std::stoull(line) == static_cast<size_t>(id), "Text block order is wrong");
        }
        std::istringstream input(serialized);
        restored.readTxtStream(input);
        check(text(restored) == serialized, "Text round trip failed");
        const Config &constConfig = restored;
        check(constConfig.ellipsoidalInclusions().size() == 2, "Const accessor failed");

        config.writeBin(0);
        std::ifstream binary(config.getBinFilename(0), std::ios::binary);
        for(auto count : counts)
        {
            size_t actual = 0;
            binary.read(reinterpret_cast<char *>(&actual), sizeof(actual));
            check(binary.good() && actual == count, "Binary count order is wrong");
        }
        model::SphericalInclusionIO<3> sphere;
        binary.read(reinterpret_cast<char *>(&sphere), sizeof(sphere));
        check(sphere.inclusionID == 10, "Binary spherical block is misplaced");
        for(size_t id : {20, 21})
        {
            model::EllipsoidalInclusionIO<3> ellipsoid;
            binary.read(reinterpret_cast<char *>(&ellipsoid), sizeof(ellipsoid));
            check(ellipsoid.inclusionID == id, "Binary ellipsoidal block is misplaced");
        }
        model::PolyhedronInclusionIO<3> polyhedron;
        binary.read(reinterpret_cast<char *>(&polyhedron), sizeof(polyhedron));
        check(polyhedron.inclusionID == 30, "Binary polyhedron block is misplaced");
        binary.close();
        restored.readBin(0);
        check(text(restored) == serialized, "Binary round trip failed");
        restored.bin2txt(0, false);
        restored.readTxt(0);
        check(text(restored) == serialized, "Binary-to-text conversion failed");

        // Trailing text matrices must remain aligned after the new inclusion block.
        config.displacementMatrix().setConstant(2, 6, 1.25);
        config.cdMatrix().setConstant(1, config.cdMatrix().cols(), 2.5);
        std::istringstream withMatrices(text(config));
        restored.readTxtStream(withMatrices);
        check(text(restored) == text(config), "Trailing matrices were shifted");

        config.clear();
        check(config.ellipsoidalInclusions().empty(), "clear() retained ellipsoids");
        const auto empty = text(config);
        std::istringstream emptyInput(empty);
        restored.readTxtStream(emptyInput);
        check(text(restored) == empty, "Empty text read retained old data");
        config.writeBin(1);
        restored.readBin(1);
        check(text(restored) == empty, "Empty binary round trip failed");
        std::filesystem::remove_all(directory);
        std::cout << "PASS DefectiveCrystal configuration defaults, inclusion ordering, text/binary round trips, conversion, matrices, and clear\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
