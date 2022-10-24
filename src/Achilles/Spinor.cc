#include "Achilles/Spinor.hh"
#include "Achilles/Utilities.hh"

using achilles::Spinor;
using achilles::SpinMatrix;

Spinor::Spinor(bool type, bool bar, const int &hel, const FourVector &mom, int ms)
        : m_type{type}, m_bar{bar}, m_hel{hel} {

    bool mode = m_type ^ (m_hel < 0);
    if(std::abs(mom.Px()) == 0.0 && std::abs(mom.Py()) == 0.0 && std::abs(mom.Pz()) == 0) {
        Complex rte(sqrt(Complex(mom.E())));
        if(mode) { // u+(p,m) or v-(p,m)
            m_u[2] = rte;
            m_u[3] = 0.0;
        } else { // u-(p,m) or v+(p,m)
            m_u[0] = 0.0;
            m_u[1] = -rte;
        }
        double sgn = m_type ? 1.0 : -1.0;
        size_t r = mode ? 0 : 2;
        m_u[0+r] = sgn*m_u[2-r];
        m_u[1+r] = sgn*m_u[3-r];
        if(m_bar) {
            m_bar = false;
            *this=Bar();
        }
    } else {
        FourVector ph(mom.E() < 0.0 ? -mom.P() : mom.P(), mom.Px(), mom.Py(), mom.Pz());
        if(mode) { // u+(p,m) or v-(p,m)
            WeylSpinor sh(true,ph);
            m_u[2] = sh.U1();
            m_u[3] = sh.U2();
        } else { // u-(p,m) or v+(p,m)
            WeylSpinor sh(false,ph);
            if(mom.E() < 0.0) sh=-sh;
            m_u[0] = sh.U2();
            m_u[1] = -sh.U1();
        }
        double m2 = mom.M2();
        if(!IsZero(m2, tol)) {
            double sgn = m_type ^ (ms < 0) ? 1.0 : -1.0;
            double omp = sqrt((mom.E() + ph.E())/(2*ph.E()));
            double omm = sqrt((mom.E() - ph.E())/(2*ph.E()));
            size_t r = mode ? 0 : 2;
            m_u[0+r] = sgn*omm*m_u[2-r];
            m_u[1+r] = sgn*omm*m_u[3-r];
            m_u[2-r] *= omp;
            m_u[3-r] *= omp;
        }
        if(m_bar) {
            m_bar = false;
            *this = Bar();
        }
    }
}

Spinor Spinor::Bar() const {
    return Spinor(m_type, !m_bar, m_hel,
                  {std::conj(m_u[2]), std::conj(m_u[3]),
                   std::conj(m_u[0]), std::conj(m_u[1])});
}

SpinMatrix Spinor::outer(const Spinor &other) const {
    if(m_bar) throw std::logic_error("Achilles::SpinorError: LHS spinor should not be barred");
    if(!other.m_bar) throw std::logic_error("Achilles::SpinorError: RHS spinor should be barred");
    if(m_type != other.m_type) throw std::logic_error("Achilles::SpinorError: Both spinors should be u or v");
    
    SpinMatrix s;

    if(m_hel != other.m_hel) return s;

    for(size_t i = 0; i < 4; ++i) {
        for(size_t j = 0; j < 4; ++j) {
            s[4*i + j] = m_u[i]*other[j];
        }
    }
    
    return s;
}

SpinMatrix SpinMatrix::PL() { return (Identity() - Gamma_5())/2.0; }
SpinMatrix SpinMatrix::PR() { return (Identity() + Gamma_5())/2.0; }

SpinMatrix SpinMatrix::Slashed(const FourVector &mom) {
    return mom.E()*Gamma_0() - mom.Px()*Gamma_1() - mom.Py()*Gamma_2() - mom.Pz()*Gamma_3();
}

SpinMatrix SpinMatrix::SigmaMuNu(size_t mu, size_t nu) {
    return std::complex<double>(0, 1)/2.0*(GammaMu(mu)*GammaMu(nu) - GammaMu(nu)*GammaMu(mu));
}

Spinor achilles::operator*(const Spinor &lhs, const SpinMatrix &rhs) {
    Spinor s(lhs.Type(), lhs.Barred(), lhs.Helicity(), std::array<std::complex<double>, 4>{0, 0, 0, 0});
    for(size_t i = 0; i < 4; ++i) {
        for(size_t j = 0; j < 4; ++j) {
            s[j] += lhs[i]*rhs[4*i+j];
        }
    }

    return s;
}

Spinor achilles::operator*(const SpinMatrix &lhs, const Spinor &rhs) {
    Spinor s(rhs.Type(), rhs.Barred(), rhs.Helicity(), std::array<std::complex<double>, 4>{0, 0, 0, 0});
    for(size_t i = 0; i < 4; ++i) {
        for(size_t j = 0; j < 4; ++j) {
            s[i] += lhs[4*i+j]*rhs[j];
        }
    }

    return s;
}
