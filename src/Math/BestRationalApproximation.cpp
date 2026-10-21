/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <BestRationalApproximation.h>
#include <math.h>

namespace model
{

    auto BestRationalApproximation::rat_approx(double f, LongIntType md) -> std::pair<LongIntType, LongIntType>
    {
        if (fabs(f) < 1.0 / md)
        {
            return std::pair<LongIntType, LongIntType>(0, 1);
        }

        /*  a: continued fraction coefficients. */
        LongIntType a, h[3] = {0, 1, 0}, k[3] = {1, 0, 0};
        LongIntType x, d, n = 1;
        long int i, neg = 0;

        if (md <= 1)
        { // *denom = 1; *num = (LongIntType) f; return;
            return std::pair<LongIntType, LongIntType>(f, 1);
        }

        if (f < 0.0)
        {
            neg = 1;
            f = -f;
        }

        while (f != floor(f))
        {
            n <<= 1;
            f *= 2;
        }

        d = f;

        /* continued fraction and check denominator each step */
        for (i = 0; i < 64; i++)
        {

            a = n ? d / n : 0;
            if (i && !a)
            {
                break;
            }

            x = d;
            d = n;

            n = x % n;

            x = a;
            if (k[1] * a + k[0] >= md)
            {
                x = (md - k[0]) / k[1];
                if (x * 2 >= a || k[1] >= md)
                {
                    i = 65;
                }
                else
                {
                    break;
                }
            }

            h[2] = x * h[1] + h[0];
            h[0] = h[1];
            h[1] = h[2];
            k[2] = x * k[1] + k[0];
            k[0] = k[1];
            k[1] = k[2];
        }
        return std::pair<LongIntType, LongIntType>(neg ? -h[1] : h[1], k[1]);
    }

    BestRationalApproximation::BestRationalApproximation(const double &f, const LongIntType &md)
        : intPair(rat_approx(f, md)), num(intPair.first), den(intPair.second)
    {
    }
} // namespace model
