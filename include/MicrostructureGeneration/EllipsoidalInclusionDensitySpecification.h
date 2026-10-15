/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_EllipsoidalInclusionDensitySpecification_H_
#define model_EllipsoidalInclusionDensitySpecification_H_
#include <MicrostructureSpecificationBase.h>
#include <Eigen/Dense>

namespace model
{
    struct EllipsoidalInclusionDensitySpecification : public MicrostructureSpecificationBase
    {
        static constexpr int dim = 3;
        double targetDensity = 0.0; // Number per cubic meter.
        Eigen::Matrix3d crystalAxes = Eigen::Matrix3d::Identity();
        // Columns of crystalAxes are orthogonal unit directions (dimensionless).
        // Arithmetic means and standard deviations of independent semiaxis lengths [m].
        Eigen::Matrix<double, 1, 3> axesLognormalDistributionMean_SI = Eigen::Matrix<double, 1, 3>::Zero();
        Eigen::Matrix<double, 1, 3> axesLognormalDistributionStd_SI = Eigen::Matrix<double, 1, 3>::Zero();
        Eigen::Matrix<double, 1, 9> crystalEigenDistortion = Eigen::Matrix<double, 1, 9>::Zero();
        Eigen::Matrix<double, 1, 3> patternVector_SI = Eigen::Matrix<double, 1, 3>::Zero();
        bool allowOverlap = false;
        bool allowOutside = false;
        double velocityReductionFactor = 1.0;
        int phaseID = 0;
        size_t maxAttempts = 10000; // Consecutive rejected placements before throwing.

        EllipsoidalInclusionDensitySpecification();
        EllipsoidalInclusionDensitySpecification(const std::string &fileName);
    };
} // namespace model
#endif
