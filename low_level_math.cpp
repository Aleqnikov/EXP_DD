#include "low_level_math.h"


double math::quick_two_sum(double a, double b, double& err)
{
    double s = a + b;
    err = b  - (s - a);

    return s;
}

double math::two_sum (double a, double b, double& err)
{
    double s = a + b;
    double v = s - a;
    err =  (a -( s - v)) + (b - v);

    return s;
}

double math::two_diff(double a, double b, double &err)
{
    double s = a  - b;
    double t = s - a;
    err = (a - (s - t)) - (b + t);
    return s;
}

void math::split(double a, double &hi, double &lo)
{
    /**
     * Нельзя просто так взять и выполнить классический алгоритм разделения, так как
     * при слишком больших значениях "a" промежуточное произведение может вызвать переполнение.
     * Поэтому добавлена проверка: если число критическое, мы сначала "сжимаем" его,
     * проводим разделение и затем возвращаем масштаб обратно.
     */

    double t;

    if (a > thresh || a < -thresh)
    {
        a *=  3.7252902984619140625e-09;    // 2^-28

        t = split_x * a;

        hi = t - (t - a);
        lo = a - hi;
        hi *= 268435456.0;                  // 2^28
        lo *= 268435456.0;                  // 2^28
    }
    else
    {
        t = split_x * a;
        hi = t - (t - a);
        lo = a - hi;
    }
}

double math::square(double a, double& err)
{
    double hi, lo;

    double q = a*a;
    split(a, hi, lo);

    err = ((hi*hi - q) + 2.0 * hi * lo) + lo * lo;

    return q;
}

double math::two_prod (double a, double b, double &err)
{
    double a_h, a_l, b_h, b_l;
    double p = a * b;
    split(a, a_h, a_l);
    split(b, b_h, b_l);

    err = ((a_h * b_h - p) + a_h * b_l + a_l * b_h) + a_l * b_l;

    return p;
}