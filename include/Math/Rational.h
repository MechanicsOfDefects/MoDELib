/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Rational_h_
#define model_Rational_h_

#include <Eigen/Core>

#include <ostream>

namespace model
{
    class Rational
    {

        /**********************************************************************/
        static int gcd(const size_t &a, const size_t &b);

        /**********************************************************************/
        static int sgn(const long int &a);

      public:
        long int n; // numerator
        long int d; // denominator

        Rational();

        Rational(const long int &n_in);

        Rational(const long int &n_in, const long int &d_in);

        double asDouble() const;

        bool operator==(const long int &other) const;

        bool operator!=(const long int &other) const;

        Rational operator*(const Rational &r2) const;

        Rational operator/(const Rational &r2) const;

        Rational operator+(const Rational &r2) const;

        Rational operator-(const Rational &r2) const;

        Rational operator*(const long int &i) const;

        Rational operator/(const long int &i) const;

        Rational operator+(const long int &i) const;

        Rational operator-(const long int &i) const;

        /*************************************************************/
        // friend T& operator <<
        // template <class T>
        friend std::ostream &operator<<(std::ostream &os, const Rational &r);
    };

    //    Rational rat(const long int& n,const long int& d)
    //    {
    //        return Rational(n,d);
    //    }

} // namespace model

namespace Eigen
{
    template <>
    struct NumTraits<model::Rational>
        : NumTraits<long int> // permits to get the epsilon, dummy_precision, lowest, highest functions
    {
        typedef long int Real;
        typedef long int NonInteger;
        typedef long int Nested;
        enum
        {
            IsComplex = 0,
            IsInteger = 1,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 3,
            MulCost = 3
        };
    };
} // namespace Eigen

#endif
