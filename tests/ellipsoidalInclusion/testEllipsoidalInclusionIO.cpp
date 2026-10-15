#include <EllipsoidalInclusionIO.h>
#include <EllipsoidalInclusion.h>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
    void check(const bool condition, const char *message)
    {
        if(!condition)
        {
            throw std::runtime_error(message);
        }
    }
} // namespace

int main()
{
    try
    {
        using IO = model::EllipsoidalInclusionIO<3>;
        using Inclusion = model::EllipsoidalInclusion<3>;
        const IO empty;
        check(empty.inclusionID == 0 && empty.C.isZero() && empty.axes.isZero() && empty.eT.isZero() &&
                  empty.mobilityReduction == 1.0 && empty.phaseID == 0,
              "Default record is not initialized");

        // Read an independently specified record to catch row/column transpositions.
        std::stringstream input("42 1 -2 3  0 -2 0  3 0 0  0 0 1 "
                                "0.01 0.003 -0.004  0.003 -0.002 0.005 "
                                "-0.004 0.005 0.007  0.625 7");
        const IO record(input);
        Eigen::Matrix3d axes;
        axes << 0, -2, 0, 3, 0, 0, 0, 0, 1;
        Eigen::Matrix3d strain;
        strain << 0.01, 0.003, -0.004, 0.003, -0.002, 0.005, -0.004, 0.005, 0.007;
        check(!input.fail() && record.inclusionID == 42 && record.phaseID == 7 && record.mobilityReduction == 0.625 &&
                  record.C.isApprox(Eigen::Vector3d(1, -2, 3)) && record.axes.isApprox(axes) &&
                  record.eT.isApprox(strain),
              "Record field order is incorrect");

        const Inclusion original(record.C, record.axes, record.eT, 0.31, 43.0, record.mobilityReduction, record.phaseID,
                                 nullptr);
        const IO fromInclusion(original);
        check(fromInclusion.inclusionID == original.sID, "Inclusion ID was not copied");
        const IO fromFields(record.inclusionID, record.C, record.axes, record.eT, record.mobilityReduction,
                            record.phaseID);
        for(const IO *source : {&fromInclusion, &fromFields})
        {
            std::stringstream text;
            text << *source;
            check(text.str().find('\n') == std::string::npos, "Output spans multiple lines");
            const IO restored(text);
            check(!text.fail() && restored.inclusionID == source->inclusionID && restored.C.isApprox(source->C) &&
                      restored.axes.isApprox(source->axes) && restored.eT.isApprox(source->eT) &&
                      restored.mobilityReduction == source->mobilityReduction && restored.phaseID == source->phaseID,
                  "Text round trip failed");
            std::string extra;
            check(!(text >> extra), "Output contains unexpected extra fields");

            const Inclusion reconstructed(restored.C, restored.axes, restored.eT, 0.31, 43.0,
                                          restored.mobilityReduction, restored.phaseID, nullptr);
            for(const Eigen::Vector3d &point : {record.C, Eigen::Vector3d(5, 6, -7)})
            {
                check(reconstructed.contains(point) == original.contains(point), "Round trip changed containment");
                check(reconstructed.stress(point).isApprox(original.stress(point), 1.e-13),
                      "Round trip changed stress");
            }
        }

        for(const char *invalid : {"", "42 1 2", "42 invalid"})
        {
            std::stringstream text(invalid);
            const IO rejected(text);
            check(text.fail(), "Invalid record did not set stream failure");
        }
        std::cout << "PASS EllipsoidalInclusionIO field order, round trips, reconstruction, and invalid input\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
