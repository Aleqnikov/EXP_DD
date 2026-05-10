#ifndef LOW_LEVEL_MATH_H
#define LOW_LEVEL_MATH_H

#include <cmath>

namespace math
{
    // Верно тогда и только тогда, когда |a| >= |b|
    inline double quick_two_sum(double a, double b, double& err)
    {
        double s = a + b;
        err = b - (s - a);
        return s;
    }

    inline double two_sum(double a, double b, double& err)
    {
        double s = a + b;
        double v = s - a;
        err = (a - (s - v)) + (b - v);
        return s;
    }

    inline double two_diff(double a, double b, double& err)
    {
        double s = a - b;
        double t = s - a;
        err = (a - (s - t)) - (b + t);
        return s;
    }


    inline double two_prod(double a, double b, double& err)
    {
        double p = a * b;
        err = std::fma(a, b, -p);
        return p;
    }

    inline double square(double a, double& err)
    {
        double q = a * a;
        err = std::fma(a, a, -q);
        return q;
    }
}

#endif //LOW_LEVEL_MATH_H