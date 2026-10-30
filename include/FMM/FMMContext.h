/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_FMMContext_H_
#define model_FMMContext_H_
#include <FMM/BoxTree.h>
#include <string>
namespace model
{
    // Shared geometry for all physics. Kernel registration will be added when
    // concrete operators are introduced; no force or energy kernels yet.
    template<int dim>
    struct FMMContext
    {
        using Vector = typename BoxTree<dim>::Vector;
        const unsigned fmmTreeDepth;
        const std::shared_ptr<BoxTree<dim>> boxTree;
        FMMContext(const Vector& lower, const Vector& upper, const std::string& inputFile);
        FMMContext(const FMMContext&) = delete;
        FMMContext& operator=(const FMMContext&) = delete;
    };
    extern template struct FMMContext<1>;
    extern template struct FMMContext<2>;
    extern template struct FMMContext<3>;
}
#endif
