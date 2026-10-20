/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <GaussLegendreRecursiveRelation.h>

namespace model
{

    auto GaussLegendreRecursiveRelation::a(const int &n) -> double
    {
        return (2.0 * n - 1.0) / n;
    }

    auto GaussLegendreRecursiveRelation::b(const int &n) -> double
    {
        return 0.0 * n;
    }

    auto GaussLegendreRecursiveRelation::c(const int &n) -> double
    {
        return (n - 1.0) / n;
    }
} // namespace model
