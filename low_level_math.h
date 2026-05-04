#ifndef LOW_LEVEL_MATH_H
#define LOW_LEVEL_MATH_H

namespace math
{
    // Важно отметить, верно тогда и только тогда когда |a| > |b|
    double quick_two_sum(double a, double b, double& err);

    double two_sum (double a, double b, double& err);

    double two_diff(double a, double b, double &err);

    // Для функции split.
    const double thresh = 6.69692879491417e+299;    // = 2^99.
                                                    // как бы верхний возможный лимит,
                                                    // после которого будет переполнение.
    const double split_x =  134217729.0;            // = 2^27 + 1,
                                                    // Как бы число которое разделяет
                                                    // число на основную часть и на его ошибку

    void split(double a, double &hi, double &lo);

    double square(double a, double& err);

    double two_prod (double a, double b, double &err);
}

#endif //LOW_LEVEL_MATH_H