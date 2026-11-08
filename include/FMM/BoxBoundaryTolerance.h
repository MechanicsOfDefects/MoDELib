#ifndef model_BoxBoundaryTolerance_H_
#define model_BoxBoundaryTolerance_H_
#include <algorithm>
#include <limits>
namespace model
{
    // Boundary projection, segment reconstruction and quadrature interpolation
    // accumulate more roundoff than a single arithmetic operation. This is a
    // numerical allowance, not a geometric extension of the root. The absolute
    // floor matches the float-epsilon distance tolerance used by DD boundary
    // geometry (coordinates are in Burgers-vector units).
    template<class Vector> double boxBoundaryTolerance(const Vector& lower,const Vector& upper)
    {
        return std::max(double(std::numeric_limits<float>::epsilon()),
            1024*std::numeric_limits<double>::epsilon()*std::max({1.0,lower.norm(),upper.norm()}));
    }
}
#endif
