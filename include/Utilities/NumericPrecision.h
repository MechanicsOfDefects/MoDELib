/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NumericPrecision_H_
#define model_NumericPrecision_H_

#include <limits>

namespace model
{
	
    template<typename T>
    struct NumericPrecision
    {
        static constexpr T epsilon=std::numeric_limits<double>::epsilon();
    };

    template<>
    struct NumericPrecision<double>
    {
        static constexpr double epsilon=std::numeric_limits<double>::epsilon();
    };

    
    template<>
    struct NumericPrecision<float>
    {
        static constexpr double epsilon=std::numeric_limits<float>::epsilon();
    };

} // namespace model
#endif
