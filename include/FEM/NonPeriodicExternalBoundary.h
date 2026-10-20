/* This file is part of MoDELib, distributed under the GNU GPL v2. */
#ifndef model_NonPeriodicExternalBoundary_H_
#define model_NonPeriodicExternalBoundary_H_
#include <ExternalBoundary.h>
namespace model
{
    // Paired faces have natural periodic traction balance, not free-surface loads.
    struct NonPeriodicExternalBoundary
    {
        template<class FE,int order,template<unsigned short,size_t> class Rule>
        static auto boundary(const FE& fe)
        {
            IntegrationDomain<FE,1,order,Rule> domain;
            for(const auto& face : ExternalBoundary::boundary<FE,order,Rule>(fe))
            {
                if(!fe.mesh.periodicSimplexPairs().count(&face.first->simplex.child(face.second)))
                {
                    domain.emplace_back(face.first,face.second);
                }
            }
            domain.bind(fe);
            return domain;
        }
    };
}
#endif
