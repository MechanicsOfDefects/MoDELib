/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearForm_H_
#define model_LinearForm_H_

#include <Eigen/Dense>
#include <EvalExpression.h>
#include <EvalFunction.h>
#include <TestExpression.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template <typename TestType,typename EvalType>
    struct LinearForm : public LinearWeakExpression<LinearForm<TestType,EvalType>>
    {

        static_assert(TestType::rows==EvalType::rows,"YOU ARE CREATING A LinearForm BETWEEN EXPRESSIONS WITH DIFFERENT NUMBER OF ROWS");

        typedef typename TestType::TrialFunctionType  TrialFunctionType;
        constexpr static int evalCols=EvalType::cols;

        typedef LinearForm<TestType,EvalType> LinearFormType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;

        constexpr static int dim=TypeTraits< TrialFunctionType>::dim;
        constexpr static int dofPerNode=TypeTraits< TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits< TrialFunctionType>::dofPerElement;

        const std::shared_ptr<TestExpression<TestType>> testExp;
        const std::shared_ptr<EvalType> evalExp;
        LinearForm(std::shared_ptr<TestExpression<TestType>> x, std::shared_ptr<EvalType> y) :
            testExp(std::move(x)), evalExp(std::move(y))
        {
            if(!testExp || !evalExp)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        Eigen::Matrix<double,dofPerElement,evalCols> operator()(const ElementType& ele,
                                                                const BaryType& bary) const
        {
            return testExp->sfm(ele,bary).transpose()*(*evalExp)(ele,bary);
        }

        const auto& field() const
        {
            return testExp->field();
        }

    };

    template<class T1,class T2, std::enable_if_t<std::is_base_of_v<EvalFunction<T2>,T2>,int> = 0>
    auto operator,(const std::shared_ptr<TestExpression<T1>>& t,const std::shared_ptr<T2>& e)
    {
        return std::make_shared<LinearForm<T1,T2>>(t,e);
    }


}
#endif
