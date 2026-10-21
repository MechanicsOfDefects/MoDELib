/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <Rational.h>

namespace model
{

    auto Rational::gcd(const size_t &a, const size_t &b) -> int
    {
        return b > 0 ? gcd(b, a % b) : a;
    }

    auto Rational::sgn(const long int &a) -> int
    {
        return a < 0 ? -1 : 1;
    }

    Rational::Rational() : n(0), d(1)
    {
    }

    Rational::Rational(const long int &n_in) : n(n_in), d(1)
    {
    }

    Rational::Rational(const long int &n_in, const long int &d_in)
        : n(sgn(d_in) * n_in / gcd(abs(n_in), abs(d_in))), d(sgn(d_in) * d_in / gcd(abs(n_in), abs(d_in)))
    {
    }

    auto Rational::asDouble() const -> double
    {
        return double(n) / double(d);
    }

    auto Rational::operator==(const long int &other) const -> bool
    {
        switch (d)
        {
        case 0:
            return false;
            break;

        case 1:
            return n == other;
            break;

        default:
            return n == 0 && other == 0;
            break;
        }
    }

    auto Rational::operator!=(const long int &other) const -> bool
    {
        return !(*this == other);
    }

    auto Rational::operator*(const Rational &r2) const -> Rational
    {
        return Rational(n * r2.n, d * r2.d);
    }

    auto Rational::operator/(const Rational &r2) const -> Rational
    {
        return Rational(n * r2.d, d * r2.n);
    }

    auto Rational::operator+(const Rational &r2) const -> Rational
    {
        return Rational(n * r2.d + r2.n * d, d * r2.d);
    }

    auto Rational::operator-(const Rational &r2) const -> Rational
    {
        return Rational(n * r2.d - r2.n * d, d * r2.d);
    }

    auto Rational::operator*(const long int &i) const -> Rational
    {
        return Rational(n * i, d);
    }

    auto Rational::operator/(const long int &i) const -> Rational
    {
        return Rational(n, d * i);
    }

    auto Rational::operator+(const long int &i) const -> Rational
    {
        return *this + Rational(i, 1);
    }

    auto Rational::operator-(const long int &i) const -> Rational
    {
        return *this - Rational(i, 1);
    }

    std::ostream &operator<<(std::ostream &os, const Rational &r)
    {
        os << r.n << "/" << r.d;
        return os;
    }
} // namespace model
