#ifndef DD_H
#define DD_H

#include <cmath>
#include "low_level_math.h"

/**
 * Содержится структура, описывающяя тип double-double.
 * а также необходимые  для работы функции.
 */

struct dd
{
    double hi;
    double lo;

    dd() : hi(0.0), lo(0.0) {}
    dd (double hi, double lo) : hi(hi), lo(lo) {};
    dd (double h) { hi = h; lo = 0.0; }

    dd operator+(const dd& other) const;
    dd operator-(const dd &other) const;
    dd operator+=(double a);
    dd operator*(const dd& other) const;
    dd operator*(const double b) const;
};

dd mul(dd a, double b);

dd square( dd a);

dd ldexp(const dd &a, int exp);

#endif //DD_H