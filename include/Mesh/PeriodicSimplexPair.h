#ifndef model_PeriodicSimplexPair_H_
#define model_PeriodicSimplexPair_H_
#include <Simplex.h>
#include <array>
#include <map>
namespace model
{
    template<int dim>
    struct PeriodicSimplexPair
    {
        const Simplex<dim,dim-1>* partner;
        Eigen::Matrix<double,dim,1> shift; // partner position = source + shift
        // Source vertex i corresponds to partner vertex vertexPermutation[i].
        std::array<size_t,dim> vertexPermutation;
    };
    template<int dim>
    using PeriodicSimplexPairMap=std::map<const Simplex<dim,dim-1>*,PeriodicSimplexPair<dim>>;
}
#endif
