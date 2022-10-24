#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <iostream>

#include "fmt/format.h"
#include "Achilles/Interpolation.hh"

constexpr double ipow(double x, size_t exponent) {
    return (exponent == 0) ? 1 :
        (exponent % 2 == 0) ? ipow(x*x, exponent/2) :
            x * ipow(x*x, (exponent-1)/2);
}

using namespace achilles;

double achilles::Polint(const std::vector<double> &x_, const std::vector<double> &y_,
                      size_t n, double x) {
    int ns = 0;
    double dift{}, dif = std::abs(x-x_[0]);
    std::vector<double> c(n);
    std::vector<double> d(n);
    for(size_t i = 0; i < n; ++i) {
        if((dift=std::abs(x-x_[i])) < dif) {
            ns = static_cast<int>(i);
            dif=dift;
        }
        c[i] = y_[i];
        d[i] = y_[i];
    }
    double y = y_[static_cast<size_t>(ns--)];
    double ho{}, hp{}, w{}, den{};
    for(size_t m = 0; m < n - 1; ++m) {
        for(size_t i = 0; i < n-m-1; ++i) {
            ho = x_[i]-x;
            hp = x_[i+m+1]-x;
            w = c[i+1]-d[i];
            if((den=ho-hp) == 0)
                throw std::runtime_error("Polint: Error in interpolation routine");
            den = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }
        y += (2*ns+1 < static_cast<int>(n-m-1) 
                ? c[static_cast<size_t>(ns+1)] : d[static_cast<size_t>(ns--)]);
    }

    return y;
}

constexpr double Interp1D::maxDeriv;

Interp1D::Interp1D(const std::vector<double> &x, const std::vector<double> &y,
                   InterpolationType mode) : kMode{mode}  {

    if(!std::is_sorted(x.begin(), x.end()))
        throw std::logic_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        throw std::logic_error("Inputs must all be unique.");
    if(x.size() != y.size())
        throw std::logic_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;
}

void Interp1D::CubicSpline(const double& derivLeft, const double& derivRight) {
    const std::size_t n = knotX.size();
    std::vector<double> u(n);
    derivs2.resize(n);

    if(derivLeft >= maxDeriv) {
        derivs2[0] = 0.0;
        u[0] = 0.0;
    } else {
        derivs2[0] = -0.5;
        u[0] = (3/(knotX[1]-knotX[0]))*((knotY[1]-knotY[0])/(knotX[1]-knotX[0])-derivLeft);
    }

    for(std::size_t i = 1; i < n-1; ++i) {
        double sig = (knotX[i]-knotX[i-1])/(knotX[i+1]-knotX[i-1]);
        double p = sig*derivs2[i-1]+2;
        derivs2[i] = (sig-1.0)/p;
        u[i] = (knotY[i+1]-knotY[i])/(knotX[i+1]-knotX[i])
            - (knotY[i]-knotY[i-1])/(knotX[i]-knotX[i-1]);
        u[i] = (6.0*u[i]/(knotX[i+1]-knotX[i-1])-sig*u[i-1])/p;
    }

    double dn{}, un{};
    if(derivRight >= maxDeriv) {
        dn = 0.0;
        un = 0.0;
    } else {
        dn = 0.5;
        un = (3/(knotX[n-1]-knotX[n-2]))*(derivRight-(knotY[n-1]-knotY[n-2])
                /(knotX[n-1]-knotX[n-2]));
    }

    derivs2[n-1]=(un-dn*u[n-2])/(dn*derivs2[n-2]+1.0);
    for(std::size_t i = n-1; i > 0; --i) {
        derivs2[i-1] = derivs2[i-1]*derivs2[i]+u[i-1];
    }

    kSplineInit = true;
}

double Interp1D::operator()(const double& x) const {
    // Ensure the interpolation is initialized first
    if(!kSplineInit && kMode == InterpolationType::CubicSpline)
        throw std::logic_error("Interpolation is not initialized!");

    // Disallow extrapolation
    if(x > knotX.back()) 
        throw std::out_of_range(fmt::format("Input ({}) greater than maximum value ({})", x, knotX.back()));
    if(x < knotX.front()) 
        throw std::out_of_range(fmt::format("Input ({}) less than minimum value ({})", x, knotX.front()));

    // Find range by binary_search
    auto idxHigh = static_cast<size_t>(std::distance(knotX.begin(), std::upper_bound(knotX.begin(), knotX.end(), x)));
    auto idxLow = idxHigh-1;

    double result = 0;
    switch(kMode) {
        case InterpolationType::NearestNeighbor:
            result = x - knotX[idxLow] < knotX[idxHigh] - x ? knotY[idxLow] : knotY[idxHigh];
            break;
        case InterpolationType::Polynomial:
            result = PolynomialInterp(x);
            break;
        case InterpolationType::CubicSpline:
            const double height = knotX[idxHigh] - knotX[idxLow];
            const double a = (knotX[idxHigh] - x)/height;
            const double b = (x - knotX[idxLow])/height;

            result = a*knotY[idxLow] + b*knotY[idxHigh] + ((ipow(a, 3) - a)*derivs2[idxLow]
                     + (ipow(b, 3) - b)*derivs2[idxHigh])*ipow(height, 2)/6.0;
            break;
    }
    return result;
}

double Interp1D::PolynomialInterp(double x) const {
    auto idx = static_cast<size_t>(std::distance(knotX.begin(), std::lower_bound(knotX.begin(), knotX.end(), x)));
    while(idx < polyOrder/2) ++idx;
    while(knotX.size() - idx < polyOrder/2+polyOrder%2) --idx;
    std::vector<double> xInterp(&knotX[idx]-polyOrder/2,
                                &knotX[idx]+polyOrder/2+polyOrder%2);
    std::vector<double> tmp(&knotY[idx]-polyOrder/2,
                            &knotY[idx]+polyOrder/2+polyOrder%2);

    return Polint(xInterp, tmp, static_cast<size_t>(polyOrder), x);
}

Interp2D::Interp2D(const std::vector<double>& x, const std::vector<double>& y,
                   const std::vector<double>& z,
                   InterpolationType mode) : kMode{mode} {
    if(!std::is_sorted(x.begin(), x.end()))
        throw std::logic_error("Inputs must be increasing.");
    if(std::adjacent_find(x.begin(), x.end()) != x.end())
        throw std::logic_error("Inputs must all be unique.");
    if(!std::is_sorted(y.begin(), y.end()))
        throw std::logic_error("Inputs must be increasing.");
    if(std::adjacent_find(y.begin(), y.end()) != y.end())
        throw std::logic_error("Inputs must all be unique.");
    if(x.size()*y.size() != z.size())
        throw std::logic_error("Input and output arrays must be the same size.");

    knotX = x;
    knotY = y;
    knotZ = z;
}

void Interp2D::BicubicSpline() {
    for(std::size_t i = 0; i < knotX.size(); ++i) {
        derivs2.emplace_back(knotY, std::vector<double>(knotZ.begin()+static_cast<int>(i*knotY.size()), knotZ.begin()+static_cast<int>((i+1)*knotY.size())));
        derivs2.back().CubicSpline();
    }

    kSplineInit = true;
}

double Interp2D::operator()(const double& x, const double& y) const {
    // Ensure the interpolation is initialized first
    if(!kSplineInit && kMode == InterpolationType::CubicSpline)
        throw std::logic_error("Interpolation is not initialized!");

    // Disallow extrapolation
    if(x > knotX.back()) 
        throw std::out_of_range(fmt::format("Input ({}) greater than maximum x value ({})", x, knotX.back()));
    if(x < knotX.front()) 
        throw std::out_of_range(fmt::format("Input ({}) less than minimum x value ({})", x, knotX.front()));
    if(y > knotY.back()) 
        throw std::out_of_range(fmt::format("Input ({}) greater than maximum y value ({})", y, knotY.back()));
    if(y < knotY.front()) 
        throw std::out_of_range(fmt::format("Input ({}) less than minimum y value ({})", y, knotY.front()));


    double result = 0;
    switch(kMode) {
        case InterpolationType::NearestNeighbor:
            result = NearestNeighbor(x, y);
            break;
        case InterpolationType::Polynomial:
            result = PolynomialInterp(x, y);
            break;
        case InterpolationType::CubicSpline:
            std::vector<double> zTmp(knotX.size()); 
            for(std::size_t i = 0; i < knotX.size(); ++i)
                zTmp[i] = derivs2[i](y);

            Interp1D interp(knotX, zTmp);
            interp.SetType(InterpolationType::CubicSpline);
            interp.CubicSpline();
            result = interp(x);
            break;
    }

    return result;
}

double Interp2D::NearestNeighbor(double x, double y) const {
    // Find range by binary_search
    auto idxHighX = static_cast<size_t>(std::distance(knotX.begin(), std::upper_bound(knotX.begin(), knotX.end(), x)));
    auto idxLowX = idxHighX-1;
    auto idxHighY = static_cast<size_t>(std::distance(knotY.begin(), std::upper_bound(knotY.begin(), knotY.end(), y)));
    auto idxLowY = idxHighY-1;
    auto idxX = x-knotX[idxLowX] < knotX[idxHighX]-x ? idxLowX : idxHighX;
    auto idxY = y-knotY[idxLowY] < knotY[idxHighY]-y ? idxLowY : idxHighY;
    return knotZ[idxY+knotY.size()*idxX];
}

double Interp2D::PolynomialInterp(double x, double y) const {
    // Find point in x direction
    auto idxX = static_cast<size_t>(std::distance(knotX.begin(), std::lower_bound(knotX.begin(), knotX.end(), x)));
    while(idxX < polyOrderX/2) ++idxX;
    while(knotX.size() - idxX < polyOrderX/2+polyOrderX%2) --idxX;
    std::vector<double> xInterp(&knotX[idxX]-polyOrderX/2,
                                &knotX[idxX]+polyOrderX/2+polyOrderX%2);

    // Find point in y direction
    auto idxY = static_cast<size_t>(std::distance(knotY.begin(), std::lower_bound(knotY.begin(), knotY.end(), y)));
    while(idxY < polyOrderY/2) ++idxY;
    while(knotY.size() - idxY < polyOrderY/2+polyOrderY%2) --idxY;
    std::vector<double> yInterp(&knotY[idxY]-polyOrderY/2,
                                &knotY[idxY]+polyOrderY/2+polyOrderY%2);

    std::vector<double> tmp(polyOrderY);
    std::vector<double> tmp2(polyOrderX);
    for(size_t i = 0; i < polyOrderX; ++i) {
        for(size_t j = 0; j < polyOrderY; ++j) {
            tmp[j] = knotZ[idxY+j-polyOrderY/2+knotY.size()*(idxX-polyOrderX/2+i)];
        }
        tmp2[i] = Polint(yInterp, tmp, polyOrderY, y);
    }
    return Polint(xInterp, tmp2, polyOrderX, x);
}
