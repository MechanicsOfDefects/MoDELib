/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <EllipsoidalInclusionDensitySpecification.h>

namespace model
{
    EllipsoidalInclusionDensitySpecification::EllipsoidalInclusionDensitySpecification()
        : MicrostructureSpecificationBase("EllipsoidalInclusion", "Density")
    {
    }

    EllipsoidalInclusionDensitySpecification::EllipsoidalInclusionDensitySpecification(const std::string &fileName)
        : MicrostructureSpecificationBase("EllipsoidalInclusion", "Density", fileName)
    {
        targetDensity = parser->readScalar<double>("targetDensity", true);
        if(targetDensity > 0.0)
        {
            crystalAxes = parser->readMatrix<double>("crystalAxes", 3, 3, true);
            axesLognormalDistributionMean_SI =
                parser->readMatrix<double>("axesLognormalDistributionMean_SI", 1, 3, true);
            axesLognormalDistributionStd_SI = parser->readMatrix<double>("axesLognormalDistributionStd_SI", 1, 3, true);
            crystalEigenDistortion = parser->readMatrix<double>("crystalEigenDistortion", 1, 9, true);
            patternVector_SI = parser->readMatrix<double>("patternVector_SI", 1, 3, true);
            allowOverlap = parser->readScalar<int>("allowOverlap", true);
            allowOutside = parser->readScalar<int>("allowOutside", true);
            velocityReductionFactor = parser->readScalar<double>("velocityReductionFactor", true);
            phaseID = parser->readScalar<int>("phaseID", true);
            maxAttempts = parser->readScalar<size_t>("maxAttempts", true);
        }
    }
} // namespace model
