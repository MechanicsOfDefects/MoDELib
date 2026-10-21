/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakDiff_H_
#define model_LinearWeakDiff_H_

#include <type_traits>
#include <LinearWeakExpression.h>
#include <memory>
#include <stdexcept>

namespace model
{

    template <typename T1,typename T2>
    struct  LinearWeakDiff : public LinearWeakExpression<LinearWeakDiff<T1,T2> >
    {

        static_assert(std::is_same<typename T1::TrialFunctionType,typename T2::TrialFunctionType>::value,"YOU ARE SUMMING TrialExpressionBase OF DIFFERENT TRIALFUNCTIONS.");
        typedef typename T1::TrialFunctionType TrialFunctionType;

        const std::shared_ptr<T1> op1;

        const std::shared_ptr<T2> op2;
        LinearWeakDiff(std::shared_ptr<T1> x, std::shared_ptr<T2> y) :
            op1(std::move(x)), op2(std::move(y))
        {
            if(!op1 || !op2)
            {
                throw std::invalid_argument("Null FEM expression operand.");
            }
        }

        Eigen::Matrix<double,Eigen::Dynamic,1> globalVector() const
        {
            if(!(op2->field().sameField(field())))
            {
                throw std::runtime_error("Cannot assemble expressions from different FEM field instances.");
            }
            return op1->globalVector()-op2->globalVector();
        }

        const auto& field() const
        {
            return op1->field();
        }

    };

    template<class T1,class T2, std::enable_if_t<std::is_base_of_v<LinearWeakExpression<T1>,T1> && std::is_base_of_v<LinearWeakExpression<T2>,T2>,int> = 0>
    auto operator-(const std::shared_ptr<T1>& x,const std::shared_ptr<T2>& y)
    {
        return std::make_shared<LinearWeakDiff<T1,T2>>(x,y);
    }

}
#endif
