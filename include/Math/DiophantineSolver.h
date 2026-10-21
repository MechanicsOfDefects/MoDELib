/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DiophantineSolver_h_
#define model_DiophantineSolver_h_

#include <Eigen/Dense>
// #include <LatticeGCD.h>
namespace model
{
    class DiophantineSolver
    {
      public:
        // Finds such x and y, that a  x + b  y = gcd(a, b)
        // https://github.com/ADJA/algos/blob/master/NumberTheory/DiophantineEquation.cpp
        static long int extended_gcd(long int a, long int b, long int &x, long int &y);

        // Solves equation a  x + b  y = c, writes answer to x and y
        static void solveDiophantine2vars(long int a, long int b, long int c, long int &x, long int &y);
        // http://mathafou.free.fr/exe_en/exedioph3.html
        // https://math.stackexchange.com/questions/514105/how-do-i-solve-a-linear-diophantine-equation-with-three-unknowns
        static void solveDiophantine3vars(const Eigen::Matrix<long int, 3, 1> &alphas, const long int t1,
                                          Eigen::Matrix<long int, 3, 1> &sols);
    };
} // namespace model

#endif
