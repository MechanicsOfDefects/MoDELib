/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BESTRATIONALAPPROXIMATION_H_
#define model_BESTRATIONALAPPROXIMATION_H_


// https://rosettacode.org/wiki/Convert_decimal_number_to_rational#C

#include <cstdint>
#include <utility>

namespace model
{

    class BestRationalApproximation
    {
      public:
        typedef int64_t LongIntType;

      private:
        static std::pair<LongIntType, LongIntType> rat_approx(double f, LongIntType md);

        const std::pair<LongIntType, LongIntType> intPair;

      public:
        const LongIntType &num;
        const LongIntType &den;

        BestRationalApproximation(const double &f, const LongIntType &md);
    };

} // namespace model
#endif
