#include "Achilles/Utilities.hh"
#include "spdlog/spdlog.h"
#include "Achilles/FourVector.hh"
#include <stdexcept>


bool achilles::CheckMasses(const std::vector<achilles::FourVector> &mom,
                           const std::vector<double> &masses, double eps) {
    if(mom.size() != masses.size()) return false;
    for(size_t i = 0; i < mom.size(); ++i) 
        if(std::abs(mom[i].M2() - masses[i]) > eps) return false;
    return true;
}

const std::array<double, 3> achilles::ToCartesian(const std::array<double, 3>& vec) {
    const double x = vec[0] * std::sin(vec[1]) * std::cos(vec[2]);
    const double y = vec[0] * std::sin(vec[1]) * std::sin(vec[2]);
    const double z = vec[0] * std::cos(vec[1]);

    return {x, y, z};
}

bool achilles::sortPairSecond(const std::pair<std::size_t, double>& a,
                              const std::pair<std::size_t, double>& b) {
    return a.second < b.second;
}

// Use the Brent algorithm to calculate the root of a given function
double achilles::Brent::CalcRoot(double a, double b) const {
    double fa = m_func(a), fb = m_func(b), fc{};
    if(std::abs(fa) < m_tol) return a;
    if(std::abs(fb) < m_tol) return b;
    if(fa*fb >= 0) throw std::domain_error("No root in given range");
    swap(fa, fb, a, b);
    double c = a;
    bool m_flag = true;
    double s{}, d = 0;
    while(fb != 0 && std::abs(b - a) > m_tol) {
        fc = m_func(c);
        if(fa != fc && fb != fc) {
            s = (a*fb*fc)/((fa-fb)*(fa-fc))+(b*fa*fc)/((fb-fa)*(fb-fc))+(c*fa*fb)/((fc-fa)*(fc-fb));
        } else {
            s = b - fb*(b-a)/(fb-fa);
        }
        if((s > std::max((3*a+b)/4, b) || s < std::min((3*a+b)/4, b)) ||
                (m_flag && std::abs(s-b) >= std::abs(b-c)/2) ||
                (!m_flag && std::abs(s-b) >= std::abs(c-d)/2) ||
                (m_flag && std::abs(b-c) < m_tol) ||
                (!m_flag && std::abs(c-d) < m_tol)) {
            s = (a+b)/2;
            m_flag = true;
        } else {
            m_flag = false;
        }

        double fs = m_func(s);
        d = c;
        c = b;
        if(fa*fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        swap(fa, fb, a, b);
    }
    return b;
}

void achilles::Brent::GetBracket(double &a, double &b, double &c,
                                 double &fa, double &fb, double &fc) const {
    fa = m_func(a);
    fb = m_func(b);

    // Ensure we are going downhill
    swap(fa, fb, a, b);

    // First guess for c
    c = b+gold*(b-a);
    fc = m_func(c);
    double u, fu;

    // Keep updating guess until we bracket, start with parabolic extrapolation
    while (fb > fc) {
        double r = (b-a)*(fb-fc);
        double q = (b-c)*(fb-fa);
        // tiny is used to prevent division by zero
        u = b-((b-c)*q-(b-a)*r)/(2.0*std::max(std::abs(q-r), tiny)*std::signbit(q-r));
        double ulim = b+glimit*(c-b);
        // Parabolic u is between b and c. Try it for the bracket
        if((b-u)*(u-c) > 0.0) {
            fu=m_func(u);
            // Minimum between b and c
            if(fu < fc) {
                a = b;
                b = u;
                fa = fb;
                fb = fu;

                return;
            // Minimum between a and u
            } else if(fu > fb) {
                c = u;
                fc = fu;

                return;
            }
            // Parabolic fit failed, use default magnifiction
            u = c + gold*(c-b);
            fu = m_func(u);
        // Parabolic fit is between c and its allowed limit
        } else if((c-u)*(u-ulim) > 0.0) {
            fu = m_func(u);
            if(fu < fc) {
                shift(b, c, u, c+gold*(c-b));
                shift(fb, fc, fu, m_func(u));
            }
        // Limit parabolic u to maximum allowed value
        } else if((u-ulim)*(ulim-c) >= 0.0) {
            u = ulim;
            fu = m_func(u);
        // Reject parabolic u, use default magnification
        } else {
            u = c+gold*(c-b);
            fu = m_func(u);
        }
        // Eliminate oldest point and continue
        shift(a, b, c, u);
        shift(fa, fb, fc, fu);
    }
}

double achilles::Brent::Minimize(const double &ax, const double &bx, const double &cx) const {
    double a = ax, b = bx, c = cx;
    double fa{}, fb{}, fc{};
    if(c > 1E98) {
        GetBracket(a, b, c, fa, fb, fc);
    }

    double d = 0.0, e = 0.0;
    double x, w, v, u;
    double fw, fv, fx, fu;

    // Initialize
    x=w=v=b;
    fw=fv=fx=fb;
    double xm = 0;
    for(size_t it = 0; it < itmax; ++it) {
        xm = 0.5*(a+b);
        double tol1 = m_tol*std::abs(x)+eps;
        double tol2 = 4.0*tol1;

        // Test for minimum
        if(std::abs(x-xm) <= (tol2-0.5*(b-a))) {
            return x;
        }

        // Construct a trial parabolic fit
        if(std::abs(e) > tol1) {
            double r = (x-w)*(fx-fv);
            double q = (x-v)*(fx-fw);
            double p = (x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if(q > 0.0) p = -p;
            q=std::abs(q);
            double etemp = e;
            e = d;

            // Determine acceptability of the parabolic fit, if fail use golden section
            if(std::abs(p) >= std::abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
                e = x >= xm ? a-x : b-x;
                d=cgold*e;
            // Use the parabolic step
            } else {
                d = p/q;
                u = x+d;
                if(u-a < tol2 || b-u < tol2) d=tol1*std::signbit(xm-x);
            }
        } else {
            e = x >= xm ? a-x : b-x;
            d = cgold*e;
        }
        u = std::abs(d) >= tol1 ? x+d : x+tol1*std::signbit(d);
        fu = m_func(u); // Only function evaluation per iteration
        if(fu <= fx) {
            if(u >= x) a=x;
            else b=x;
            shift(v, w, x, u);
            shift(fv, fw, fx, fu);
        } else {
            if(u < x) a=u;
            else b=u;
            if(fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if(fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }

    spdlog::warn("Brent failed to converge to minimum: {}, {}", x, std::abs(x-xm));
    return x;
}
