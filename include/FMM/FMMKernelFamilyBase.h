#ifndef model_FMMKernelFamilyBase_H_
#define model_FMMKernelFamilyBase_H_
#include <utility>
namespace model
{
    template<class Derived> class FMMKernelFamilyBase
    {
    public:
        const bool isSplit;
        explicit FMMKernelFamilyBase(bool split):isSplit(split) {}
        template<class Function> decltype(auto) withFarKernel(Function&& function) const
        {
            const auto& derived=static_cast<const Derived&>(*this);
            if(isSplit)
            {
                return std::forward<Function>(function)(derived.longKernel());
            }
            return std::forward<Function>(function)(derived.fullKernel());
        }
        template<class Function> decltype(auto) withNearKernel(Function&& function) const
        {
            const auto& derived=static_cast<const Derived&>(*this);
            if(isSplit)
            {
                return std::forward<Function>(function)(derived.shortKernel());
            }
            return std::forward<Function>(function)(derived.fullKernel());
        }
    };
}
#endif
