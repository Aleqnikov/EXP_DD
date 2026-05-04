#include <iostream>
#include <limits>
#include <cmath>


static const double inf = std::numeric_limits<double>::infinity();

// Простая функция для быстрой суммы двух чисел с учетом ошибки. |a| > |b|
double quick_two_sum(double a, double b, double& err)
{
    double s = a + b;
    err = b  - (s - a);

    return s;
}

// Простая функция для суммы двух чисел с учетом ошибки.
double two_sum (double a, double b, double& err)
{
    double s = a + b;
    double v = s - a;
    err =  (a -( s - v)) + (b - v);

    return s;
}

double two_diff(double a, double b, double &err)
{
    double s = a  - b;
    double t = s - a;
    err = (a - (s - t)) - (b + t);
    return s;
}

const double thresh = 6.69692879491417e+299; // = 2^99
const double split_x =  134217729.0;            // = 2^27 + 1

void split(double a, double &hi, double &lo)
{
    /**
     * Неправильно выполнить классический алгоритм разделения, так как
     * число a может быть слишком большим, что может вызвать переполнение,
     * поэтому, добавляется проверка, но тоЮ чтобы число не вызвало переполнение.
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

double square(double a, double& err)
{
    double hi, lo;

    double q = a*a;
    split(a, hi, lo);

    err = ((hi*hi - q) + 2.0 * hi * lo) + lo * lo;

    return q;
}

double two_prod (double a, double b, double &err)
{
    double a_h, a_l, b_h, b_l;
    double p = a * b;
    split(a, a_h, a_l);
    split(b, b_h, b_l);

    err = ((a_h * b_h - p) + a_h * b_l + a_l * b_h) + a_l * b_l;

    return p;
}



struct dd_
{
    double hi;
    double lo;

    dd_() : hi(0.0), lo(0.0) {}
    dd_ (double hi, double lo) : hi(hi), lo(lo) {};
    dd_ (double h) { hi = h; lo = 0.0; }

    dd_ operator+(const dd_& other) const
    {
        /**
         * Так называемый стиль IEEE
         */

        double s1, s2, t1, t2;

        s1 = two_sum(this->hi, other.hi, s2);
        t1 = two_sum(this->lo, other.lo, t2);

        s2 += t1;
        s1 = quick_two_sum(s1, s2, s2);

        s2 += t1;
        s1 = quick_two_sum(s1, s2, s2);

        return dd_(s1, s2);
    }

    dd_ operator-(const dd_ &other) const
    {
        double s1, s2, t1, t2;

        s1 = two_diff(this->hi, other.hi, s2);
        t1 = two_diff(this->lo, other.lo, t2);

        s2 += t1;

        s1 = quick_two_sum(s1, s2, s2);
        s2 += t2;

        s1 = quick_two_sum(s1, s2, s2);


        return dd_(s1, s2);

    }


    dd_ operator+=(double a)
    {
        double s1, s2;
        s1 = two_sum(this->hi, a, s2);
        s2 += this->lo;
        this->hi = quick_two_sum(s1, s2, this->lo);
        return *this;
    }

    dd_ operator*(const dd_& other) const
    {
        double p1, p2;

        p1 = two_prod(this->hi, other.hi, p2);
        p2 += (this->hi * other.lo + this->lo * other.hi);
        p1 = quick_two_sum(p1, p2, p2);
        return dd_(p1, p2);
    }

    dd_ operator*(const double b) const
    {
        double p1, p2;

        p1 = two_prod(this->hi, b, p2);
        p2 += (this->lo * b);
        p1 = quick_two_sum(p1, p2, p2);
        return dd_(p1, p2);
    }
};

dd_ mul_pwr2(dd_ a, double b)
{
    return dd_(a.hi * b, a.lo * b);
}

dd_ squre_two( dd_ a)
{
    double p1, p2;
    double s1, s2;

    p1 = square(a.hi, p2);
    p2 +=  2.0 * a.hi * a.lo;
    p2 += a.lo * a.lo;

    s1 = quick_two_sum(p1, p2, s2);

    return dd_(s1, s2);
}

dd_ ldexp(const dd_ &a, int exp) {
  return dd_(std::ldexp(a.hi, exp), std::ldexp(a.lo, exp));
}

namespace Consts
{
    const dd_ EXP(2.718281828459045091e+00, 1.445646891729250158e-16);
    const double  EPS = 4.93038065763132e-32;
    const dd_ PI(3.1415926535, 8.979323846e-11);
    const dd_ INF(inf, inf);
    const dd_ LN2(6.931471805599452862e-01, 2.319046813846299558e-17);
    const dd_ factorial[15] = {
        dd_ (1.66666666666666657e-01,  9.25185853854297066e-18),
        dd_ { 4.16666666666666644e-02,  2.31296463463574266e-18},
        dd_ { 8.33333333333333322e-03,  1.15648231731787138e-19},
        dd_ { 1.38888888888888894e-03, -5.30054395437357706e-20},
        dd_ { 1.98412698412698413e-04,  1.72095582934207053e-22},
        dd_ { 2.48015873015873016e-05,  2.15119478667758816e-23},
        dd_ { 2.75573192239858925e-06, -1.85839327404647208e-22},
        dd_ { 2.75573192239858883e-07,  2.37677146222502973e-23},
        dd_ { 2.50521083854417202e-08, -1.44881407093591197e-24},
        dd_ { 2.08767569878681002e-09, -1.20734505911325997e-25},
        dd_ { 1.60590438368216133e-10,  1.25852945887520981e-26},
        dd_ { 1.14707455977297245e-11,  2.06555127528307454e-28},
        dd_ { 7.64716373181981641e-13,  7.03872877733453001e-30},
        dd_ { 4.77947733238738525e-14,  4.39920548583408126e-31},
        dd_ { 2.81145725434552060e-15,  1.65088427308614326e-31}
    };

}

dd_ exp( const dd_ &a)
{
    /**
     * Редуцируем экспоненту, для того чтобы ряд тейлора быстрее сходился,
     * и затем считаем всю функцию
     */

    double k = 512.0;

    double k_inv = 1.0 / k;

    if (a.hi <= -709.0) return 0.0;
    if (a.hi >= 709.0) return Consts::INF;

    if (a.hi == 0) return 1.0;
    if (a.hi == 1.0 && a.lo == 0.0) return Consts::EXP;

    /**
     * EXP(a) = EXP(kr + m * ln(2)) = 2^m * EXP(r)^k
     *
     * Нужно подобрать таким образом m, чтобы |kr| < ln(2)/2 = 0.347
     */

    double m = std::floor(a.hi / Consts::LN2.hi + 0.5);    // + 0.5 необходим для правильного округления.
    dd_ r = mul_pwr2(a -  Consts::LN2 * m , k_inv);

    dd_ s, t, p;

    p = squre_two(r);

    s = r + mul_pwr2(p, 0.5);
    p = p * r;

    t = p * Consts::factorial[0];

    int i = 0;
    do {
        s = s + t;
        p = p * r;

        i++;

        t = p * Consts::factorial[i];
    } while (std::abs(t.hi) > k_inv * Consts::EPS && i < 15);

    s = s + t;

    // Возводим в степень.

    /**
     * Так как мы наши s = e^r - 1 то теперь мы хотим про помощи алгоритма быстрого
     * возведения в степень получить e(r)^512 то (1 + s)^2  - 1 = 1 + 2s + s^2 = 2s + s^2
     */

    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);
    s = mul_pwr2(s, 2.0) + squre_two(s);

    s += 1.0;

    return ldexp(s, static_cast<int>(m));
}

int main() {

    dd_ x(154.3456456, 0.0);
    dd_ res = exp(x);

    // Устанавливаем точность вывода
    std::cout.precision(17);
    std::cout << "e^1 = " << res.hi << " (hi)" << std::endl;
    std::cout << "err = " << res.lo << " (lo)" << std::endl;

    // Сравним с библиотечным значением
    std::cout << "std::exp(1) = " << std::exp(154.3456456) << std::endl;

    return 0;
}