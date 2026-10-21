/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TrialExpressionBase_H_
#define model_TrialExpressionBase_H_

#include <memory>
#include <type_traits>
#include <iostream>
#include <iomanip>
#include <TypeTraits.h>
#include <TestExpression.h>
#include <EvalExpression.h>

namespace model
{

    template<typename Derived>
    struct TrialExpressionBase
    {
        TrialExpressionBase()= default;
        TrialExpressionBase(const TrialExpressionBase<Derived>&) = delete;
        TrialExpressionBase(TrialExpressionBase<Derived>&&) = delete;

        Derived& derived()
        {
            return *static_cast<Derived*>(this);
        }

        Derived* p_derived()
        {
            return  static_cast<Derived*>(this);
        }

        const Derived& derived() const
        {
            return *static_cast<const Derived*>(this);
        }

        const Derived* p_derived() const
        {
            return  static_cast<const Derived*>(this);
        }

        operator const Derived&() const
        {
            return derived();
        }

    };


    template<class T, std::enable_if_t<std::is_base_of_v<TrialExpressionBase<T>,T>,int> = 0>
    auto test(const std::shared_ptr<T>& t)
    {
        return std::make_shared<TestExpression<T>>(t);
    }
    template<class T, std::enable_if_t<std::is_base_of_v<TrialExpressionBase<T>,T>,int> = 0>
    auto eval(const std::shared_ptr<T>& t)
    {
        return std::make_shared<EvalExpression<T>>(t);
    }
}
#endif
