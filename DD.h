#ifndef LB2_DD_H
#define LB2_DD_H


#include <cmath>
/**
 * Содержится структура, описывающяя тип double-double.
 * а также необходимые функции, для работы.
 */

#include "low_level_math.h"

struct dd_
{
    double hi;
    double lo;

    dd_() : hi(0.0), lo(0.0) {}
    dd_ (double hi, double lo) : hi(hi), lo(lo) {};
    dd_ (double h) { hi = h; lo = 0.0; }

    dd_ operator+(const dd_& other) const;
    dd_ operator-(const dd_ &other) const;
    dd_ operator+=(double a);
    dd_ operator*(const dd_& other) const;
    dd_ operator*(const double b) const;

};

dd_ mul_pwr2(dd_ a, double b);

dd_ squre_two( dd_ a);

dd_ ldexp(const dd_ &a, int exp);



#endif //LB2_DD_H