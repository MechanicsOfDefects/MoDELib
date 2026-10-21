/* This file is part of MoDELib, distributed under the GNU GPL v2. */
#ifndef model_TrialSum_H_
#define model_TrialSum_H_
#include <TrialExpressionBase.h>
#include <memory>
#include <stdexcept>
#include <type_traits>
namespace model
{
    template<class T1,class T2>
    struct TrialSum : TrialExpressionBase<TrialSum<T1,T2>>
    {
        using TrialFunctionType=typename T1::TrialFunctionType;
        using ElementType=typename T1::ElementType;
        using BaryType=typename T1::BaryType;
        using ShapeFunctionMatrixType=typename T1::ShapeFunctionMatrixType;
        static constexpr int rows=T1::rows;
        static constexpr int cols=T1::cols;
        static_assert(std::is_same_v<TrialFunctionType,typename T2::TrialFunctionType>);
        static_assert(rows==T2::rows && cols==T2::cols);
        const std::shared_ptr<T1> op1;
        const std::shared_ptr<T2> op2;

        TrialSum(std::shared_ptr<T1> x,std::shared_ptr<T2> y) : op1(std::move(x)),op2(std::move(y))
        {
            if(!op1 || !op2)
            {
                throw std::invalid_argument("Null FEM sum operand.");
            }
            if(!op1->field().sameField(op2->field()))
            {
                throw std::invalid_argument("Cannot add trial expressions of different FEM fields.");
            }
        }
        const auto& field() const
        {
            return op1->field();
        }
        ShapeFunctionMatrixType sfm(const ElementType& e,const BaryType& b) const
        {
            return op1->sfm(e,b)+op2->sfm(e,b);
        }
        auto sfmGrad(const ElementType& e,const BaryType& b) const
        {
            return (op1->sfmGrad(e,b)+op2->sfmGrad(e,b)).eval();
        }
        auto sfmDef(const ElementType& e,const BaryType& b) const
        {
            return (op1->sfmDef(e,b)+op2->sfmDef(e,b)).eval();
        }
    };
    template<class T1,class T2, std::enable_if_t<std::is_base_of_v<TrialExpressionBase<T1>,T1> && std::is_base_of_v<TrialExpressionBase<T2>,T2>,int> = 0>
    auto operator+(const std::shared_ptr<T1>& x,const std::shared_ptr<T2>& y)
    {
        return std::make_shared<TrialSum<T1,T2>>(x,y);
    }
}
#endif
