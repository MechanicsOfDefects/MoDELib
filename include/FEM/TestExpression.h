/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TestExpression_H_
#define model_TestExpression_H_

#include <utility>
#include <memory>
#include <stdexcept>

namespace model
{

    template<typename TrialExpressionType>
    struct TestExpression
    {

        typedef TestExpression<TrialExpressionType> TestExpressionType;
        typedef typename TrialExpressionType::ShapeFunctionMatrixType     ShapeFunctionMatrixType;
        typedef typename TrialExpressionType::BaryType BaryType;
        typedef typename TrialExpressionType::ElementType ElementType;

        const std::shared_ptr<TrialExpressionType> trialExp;

        explicit TestExpression(std::shared_ptr<TrialExpressionType> trialExp_in) : trialExp(std::move(trialExp_in))
        {
            if(!trialExp)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        ShapeFunctionMatrixType sfm(const ElementType& ele,
                                           const BaryType& bary) const
        {
            return trialExp->sfm(ele,bary);
        }

        const auto& field() const
        {
            return trialExp->field();
        }

    };

}
#endif
