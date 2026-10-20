/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearForm_H_
#define model_BilinearForm_H_

#include <type_traits>
#include <TrialExpressionBase.h>

#include <TerminalColors.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template <typename T1,typename T2>
	struct BilinearForm
    {

        static_assert(std::is_same<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE CREATING A BilinearForm OF DIFFERENT TrialFunction TYPES.");
        static_assert(T1::rows==T2::rows,"YOU ARE CREATING A BilinearForm BETWEEN A TrialExpression AND A TestExpression WITH DIFFERENT NUMBER OF ROWS");

        typedef typename TestExpression<T1>::TestExpressionType TestExpressionType;
        typedef T2 TrialExpressionType;

        typedef typename T2::TrialFunctionType TrialFunctionType;
        typedef typename TypeTraits<TrialFunctionType>::ElementType ElementType;
        typedef typename TypeTraits<TrialFunctionType>::BaryType BaryType;
        typedef typename TypeTraits<TrialFunctionType>::FiniteElementType FiniteElementType;

        constexpr static int dim=TypeTraits<TrialFunctionType>::dim;
        constexpr static int nodesPerElement=TypeTraits<TrialFunctionType>::nodesPerElement;
        constexpr static int dofPerNode=TypeTraits<TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits<TrialFunctionType>::dofPerElement;

        typedef Eigen::Matrix<double,dofPerElement,dofPerElement> ElementMatrixType;

        const std::shared_ptr<TestExpression<T1>> testExpr;
        const std::shared_ptr<TrialExpressionType> trialExpr;
        BilinearForm(std::shared_ptr<TestExpression<T1>> x, std::shared_ptr<T2> y) :
            testExpr(std::move(x)), trialExpr(std::move(y))
        {
            if(!testExpr || !trialExpr)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        ElementMatrixType operator()(const ElementType& ele, const BaryType& bary) const
        {

            if(!(testExpr->field().sameField(field())))
            {
                throw std::runtime_error("Cannot assemble expressions from different FEM field instances.");
            }
            return testExpr->sfm(ele,bary).transpose()*trialExpr->sfm(ele,bary);
        }

        const auto& field() const
        {
            return trialExpr->field();
        }

    };

    template<class T1,class T2, std::enable_if_t<std::is_base_of_v<TrialExpressionBase<T2>,T2>,int> = 0>
    auto operator,(const std::shared_ptr<TestExpression<T1>>& t,const std::shared_ptr<T2>& e)
    {
        return std::make_shared<BilinearForm<T1,T2>>(t,e);
    }

}
#endif
