#ifndef UNITS_HH
#define UNITS_HH

#include <cmath>
namespace achilles {
    constexpr double base_to_femto = 1.0e15;

    // Time literals
    constexpr double operator"" _s(long double x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _s(unsigned long long int x) {
        return static_cast<double>(x);
    }

    // Distance literals
    constexpr double operator"" _fm(long double x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _fm(unsigned long long int x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _m(long double x) {
        return static_cast<double>(x)*base_to_femto;
    }
    constexpr double operator"" _m(unsigned long long int x) {
        return static_cast<double>(x)*base_to_femto;
    }

    // Cross-section literals
    constexpr double operator"" _mb(long double x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _mb(unsigned long long int x) {
        return static_cast<double>(x);
    }

    // Energy literals
    constexpr double operator"" _MeV(long double x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _MeV(unsigned long long int x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _GeV(long double x) {
        return static_cast<double>(x*1000);
    }
    constexpr double operator"" _GeV(unsigned long long int x) {
        return static_cast<double>(x)*1000;
    }
    constexpr double operator"" _J(long double x) {
        return static_cast<double>(x)/1.602176634e-13;
    }

    // Angle literals
    constexpr double operator"" _rad(long double x) {
        return static_cast<double>(x);
    }
    constexpr double operator"" _deg(long double x) {
        constexpr double ToRads = M_PI/180;
        return static_cast<double>(x)*ToRads;
    }
}

#endif
