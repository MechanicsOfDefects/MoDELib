/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialProd_H_
#define model_TrialProd_H_

#include <utility>
#include <TypeTraits.h>
#include <Simplex.h>
#include <TrialBase.h>
#include <TrialGrad.h>
#include <TrialExpressionBase.h>
#include <Constant.h>

#include <EvalFunction.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template<int rows1,int cols1, int rows2>
    struct TrialProdRows
    {
        static_assert(cols1==rows2,"YOU ARE MULTIPLYING TrialExpressionBaseS WITH DIFFERENT SIZE.");
        constexpr static int rows=rows1;
    };

    template<int rows2>
    struct TrialProdRows<1,1,rows2>
    {
        constexpr static int rows=rows2;
    };

    template<>
    struct TrialProdRows<1,1,1>
    {
        constexpr static int rows=1;
    };

    template <typename T1, typename T2>
    struct TrialProd :
     public TrialExpressionBase<TrialProd<T1,T2> >
    {
        typedef typename T2::TrialFunctionType TrialFunctionType;
        typedef TrialBase<TrialFunctionType> TrialBaseType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;
        constexpr static int rows=TrialProdRows<T1::rows,T1::cols,T2::rows>::rows;
        constexpr static int cols=T2::cols;

        typedef Eigen::Matrix<double,rows,dofPerElement> ShapeFunctionMatrixType;
        typedef Eigen::Matrix<double,rows*dim,dofPerElement> ShapeFunctionGradMatrixType;
        typedef Eigen::Matrix<double,rows,1> EvalMatrixType;

        const std::shared_ptr<T1> op1;
        const std::shared_ptr<T2> op2;
        TrialProd(std::shared_ptr<T1> x, std::shared_ptr<T2> y) :
            op1(std::move(x)), op2(std::move(y))
        {
            if(!op1 || !op2)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        ShapeFunctionMatrixType sfm(const ElementType& ele,
                                    const BaryType& bary) const
        {
            return (*op1)(ele,bary)*op2->sfm(ele,bary);
        }

        ShapeFunctionGradMatrixType sfmGrad(const ElementType& ele,
                                                                       const BaryType& bary) const
        {
            const auto t1((*op1)(ele,bary));
            if constexpr(T1::rows==1 && T1::cols==1)
            {
                if constexpr(std::is_arithmetic_v<decltype(t1)>)
                {
                    return t1*op2->sfmGrad(ele,bary);
                }
                else
                {
                    return t1(0,0)*op2->sfmGrad(ele,bary);
                }
            }
            else
            {
                Eigen::Matrix<double,T1::rows*dim,TrialGrad<T2>::rows> t1Ex(Eigen::Matrix<double,T1::rows*dim,TrialGrad<T2>::rows>::Zero());
                for(int I=0;I<T1::rows;++I)
                {
                    for(int J=0;J<T1::cols;++J)
                    {
                        for(int d=0;d<dim;++d)
                        {
                            const int i(I*dim+d);
                            const int j(J*dim+d);
                            t1Ex(i,j)=t1(I,J);
                        }
                    }
                }
                return t1Ex*op2->sfmGrad(ele,bary);
            }
        }
        const auto& field() const
        {
            return op2->field();
        }

    };

    template<class E,class T, std::enable_if_t<std::is_base_of_v<EvalFunction<E>,E> && std::is_base_of_v<TrialExpressionBase<T>,T>,int> = 0>
    auto operator*(const std::shared_ptr<E>& e,const std::shared_ptr<T>& t)
    {
        return std::make_shared<TrialProd<E,T>>(e,t);
    }
    template<class T>
    auto operator*(double c,const std::shared_ptr<T>& t)
    {
        return make_constant(c)*t;
    }
    template<int R,int C,class T>
    auto operator*(const Eigen::Matrix<double,R,C>& c,const std::shared_ptr<T>& t)
    {
        return make_constant(c)*t;
    }

}
#endif
