/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_FMMContext_H_
#define model_FMMContext_H_
#include <FMM/BoxTree.h>
#include <FMM/FMMManager.h>
#include <string>
namespace model
{
    enum class MultipoleMode { none, ewaldFMM, FMM };

    // Shared geometry and host registry for independently registered physics kernels.
    template<int dim>
    struct FMMContext
    {
        using Vector = typename BoxTree<dim>::Vector;
        const std::shared_ptr<FMMManager> manager=std::make_shared<FMMManager>();
        const MultipoleMode multipoleMode;
        bool usesFMM() const { return multipoleMode != MultipoleMode::none; }
        bool usesEwaldFMM() const { return multipoleMode == MultipoleMode::ewaldFMM; }
        const int interpolationOrder;
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
