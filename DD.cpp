#include "DD.h"

dd dd::operator+(const dd& other) const
{
    double s1, s2, t1, t2;

    s1 = math::two_sum(this->hi, other.hi, s2);
    t1 = math::two_sum(this->lo, other.lo, t2);

    s2 += t1;
    s1 = math::quick_two_sum(s1, s2, s2);

    s2 += t2;
    s1 = math::quick_two_sum(s1, s2, s2);

    return dd(s1, s2);
}

dd dd::operator-(const dd &other) const
{
    double s1, s2, t1, t2;

    s1 = math::two_diff(this->hi, other.hi, s2);
    t1 = math::two_diff(this->lo, other.lo, t2);

    s2 += t1;

    s1 = math::quick_two_sum(s1, s2, s2);
    s2 += t2;

    s1 = math::quick_two_sum(s1, s2, s2);

    return dd(s1, s2);
}


dd dd::operator+=(double a)
{
    double s1, s2;

    s1 = math::two_sum(this->hi, a, s2);
    s2 += this->lo;
    this->hi = math::quick_two_sum(s1, s2, this->lo);
    return *this;
}

dd dd::operator*(const dd& other) const
{
    double p1, p2;

    p1 = math::two_prod(this->hi, other.hi, p2);
    p2 += (this->hi * other.lo + this->lo * other.hi);
    p1 = math::quick_two_sum(p1, p2, p2);
    return dd(p1, p2);
}

dd dd::operator*(const double b) const
{
    double p1, p2;

    p1 = math::two_prod(this->hi, b, p2);
    p2 += (this->lo * b);
    p1 = math::quick_two_sum(p1, p2, p2);
    return dd(p1, p2);
}


dd mul(dd a, double b)
{
    return dd(a.hi * b, a.lo * b);
}

dd square( dd a)
{
    double p1, p2;
    double s1, s2;

    p1 = math::square(a.hi, p2);
    p2 +=  2.0 * a.hi * a.lo;
    p2 += a.lo * a.lo;

    s1 = math::quick_two_sum(p1, p2, s2);

    return dd(s1, s2);
}

dd ldexp(const dd &a, int exp)
{
    return dd(std::ldexp(a.hi, exp), std::ldexp(a.lo, exp));
}
