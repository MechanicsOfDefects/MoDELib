/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BilinearWeakExpression_H_
#define model_BilinearWeakExpression_H_

namespace model
{

    template <typename Derived>
	struct BilinearWeakExpression
    {

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

}
#endif
