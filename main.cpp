#include <iostream>
#include <limits>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

#include "DD.h"

static const double inf = std::numeric_limits<double>::infinity();

namespace Consts
{
    const dd EXP(2.718281828459045091e+00, 1.445646891729250158e-16);
    const double EPS = 4.93038065763132e-32;    // Машинное эпсилон для dd: eps_dd = eps_d^2, eps_d = 2^{-52}
    const dd INF(inf, inf);
    const dd LN2(6.931471805599452862e-01, 2.319046813846299558e-17);

    // Коэффициенты ряда Тейлора: 1/n! для n = 3-17
    const dd factorial[15] = {
        dd (1.66666666666666657e-01,  9.25185853854297066e-18),
        dd { 4.16666666666666644e-02,  2.31296463463574266e-18},
        dd { 8.33333333333333322e-03,  1.15648231731787138e-19},
        dd { 1.38888888888888894e-03, -5.30054395437357706e-20},
        dd { 1.98412698412698413e-04,  1.72095582934207053e-22},
        dd { 2.48015873015873016e-05,  2.15119478667758816e-23},
        dd { 2.75573192239858925e-06, -1.85839327404647208e-22},
        dd { 2.75573192239858883e-07,  2.37677146222502973e-23},
        dd { 2.50521083854417202e-08, -1.44881407093591197e-24},
        dd { 2.08767569878681002e-09, -1.20734505911325997e-25},
        dd { 1.60590438368216133e-10,  1.25852945887520981e-26},
        dd { 1.14707455977297245e-11,  2.06555127528307454e-28},
        dd { 7.64716373181981641e-13,  7.03872877733453001e-30},
        dd { 4.77947733238738525e-14,  4.39920548583408126e-31},
        dd { 2.81145725434552060e-15,  1.65088427308614326e-31}
    };
}

/**
 * Быстрое возведение в степень для выражения (1 + s)^power - 1
 * с использованием формулы (1 + s)^2 - 1 = 2s + s^2.
 *
 * ВАЖНО: параметр power обязан быть степенью двойки (например, 512).
 */
dd fast_taylor_pow(dd s, int power)
{
    while (power > 1) {
        s = mul(s, 2.0) + square(s);
        power >>= 1;
    }
    return s;
}

/**
 * Вычисление exp(a) в формате double-double.
 *
 * Алгоритм:
 *   1. Аргумент преобразуется по формуле: a = m * ln(2) + k * r,
 *      откуда r = (a - m * ln(2)) / k.
 *   2. exp(r) - 1 вычисляется рядом Тейлора, важно, что быстро сходится при малом r.
 *   3. Результат восстанавливается: exp(a) = 2^m * exp(r)^k.
 */
dd exp_dd(const dd& a)
{
    const int k = 512;                  // Обязательно k = 2^p

    /**
     * Граница аргумента: максимальное значение double равно ~2^1024,
     * а его натуральный логарифм примерно равен 709.
     * При |a| > 709 результат уходит в переполнение - inf или в ноль.
     */
    const double limit_exp_arg = 709.0;

    if (a.hi <= -limit_exp_arg) return 0.0;
    if (a.hi >= limit_exp_arg)  return Consts::INF;

    if (a.hi == 0) return 1.0;
    if (a.hi == 1.0 && a.lo == 0.0) return Consts::EXP;

    /**
     * Выбор k = 512 = 2^9 - это некая константа, что сохраняет точность, и скорость быстродействия алгоритма:
     *     После вычитания m * ln(2) остаток лежит в [-ln(2)/2, ln(2)/2].
     *     Деление на k = 512 даёт столь малое r, что ряд Тейлора сходится
     *     до нужной точности примерно за 15 итераций.
     *     Степень двойки позволяет вычислить exp(r)^512 последовательными
     *     возведениями в квадрат за 9 шагов.
     *     k = 256: ряд сходился бы медленнее; k = 1024: ошибка накапливалась
     *     бы сильнее из-за слишком малого r.
     */

    double k_inv = 1.0 / k;

    /**
     * EXP(a) = EXP(k*r + m * ln(2)) = 2^m * EXP(r)^k
     */

    double m = std::floor(a.hi / Consts::LN2.hi + 0.5);    // + 0.5 необходим для правильного округления.
    dd r = mul(a - Consts::LN2 * m, k_inv);

    // Ряд Тейлора для exp(r) - 1, начиная со слагаемого r^2/2!:
    //   s = r + r^2/2! + r^3/3! + ...
    dd s, t, p;

    p = square(r);          		// p = r^2
    s = r + mul(p, 0.5);    		// s = r + r^2/2!
    p = p * r;              		// p = r^3

    t = p * Consts::factorial[0];   // t = r^3 / 3!

    int i = 0;
    do {
        s = s + t;
        p = p * r;

        i++;

        t = p * Consts::factorial[i];
    } while (std::abs(t.hi) > k_inv * Consts::EPS && i < 15 - 1);

    s = s + t;

    // Возводим в степень.

    /**
     * Теперь s = exp(r) - 1.
     * Применяем быстрое возведение в степень:
     * (1 + s)^k - 1, используя (1 + s)^2 - 1 = 2s + s^2.
     * После этого прибавляем 1 и получаем exp(r)^k = exp(k*r).
     */

    s = fast_taylor_pow(s, k);

    s += 1.0;

    return ldexp(s, static_cast<int>(m));
}


void print_dd(const dd& x)
{
    if (std::isinf(x.hi) || std::isnan(x.hi)) {
        std::cout << x.hi << std::endl;
        return;
    }

    const int precision = 32;

    int exp_hi = 0;
    if (x.hi != 0.0) {
        exp_hi = static_cast<int>(std::floor(std::log10(std::abs(x.hi))));
    }


    double lo_scaled = x.lo * std::pow(10.0, -exp_hi);

    std::ios_base::fmtflags fl = std::cout.flags();

    std::cout << std::scientific << std::setprecision(precision);
    std::cout << x.hi  << std::endl;


    std::cout << std::fixed << std::setprecision(precision);
    std::cout << lo_scaled << "e" << (exp_hi >= 0 ? "+" : "") << exp_hi  << std::endl;

    std::cout.flags(fl);
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " -x <hi_value> [-l <lo_value>]" << std::endl;
        return 1;
    }

    double x_hi = 0.0;
    double x_lo = 0.0;
    bool hi_provided = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-x" || arg == "--xi") && i + 1 < argc) {
            x_hi = std::stod(argv[++i]);
            hi_provided = true;
        } else if ((arg == "-l" || arg == "--lo") && i + 1 < argc) {
            x_lo = std::stod(argv[++i]);
        }
    }

    if (!hi_provided) {
        std::cerr << "Error: Argument -x is required." << std::endl;
        return 1;
    }

    dd x(x_hi, x_lo);
    dd res = exp_dd(x);

    print_dd(res);

    return 0;
}