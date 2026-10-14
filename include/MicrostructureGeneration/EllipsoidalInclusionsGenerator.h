/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_EllipsoidalInclusionsGenerator_H_
#define model_EllipsoidalInclusionsGenerator_H_
#include <EllipsoidalInclusionIndividualSpecification.h>
#include <EllipsoidalInclusionDensitySpecification.h>
namespace model
{
    class MicrostructureGenerator;
    class EllipsoidalInclusionsGenerator
    {
    public:
        EllipsoidalInclusionsGenerator(const EllipsoidalInclusionDensitySpecification& spec,
                                      MicrostructureGenerator& mg);
        EllipsoidalInclusionsGenerator(const EllipsoidalInclusionIndividualSpecification& spec,
                                      MicrostructureGenerator& mg);
    };
}
#endif
