#ifndef SPINOR_HH
#define SPINOR_HH

// Implementation of spinor classes closely related to that in Sherpa

#include <complex>
#include <iostream>

#include "Achilles/FourVector.hh"

namespace achilles {

class WeylSpinor {
    public:
        using Complex = std::complex<double>;

        WeylSpinor() = default;
        WeylSpinor(bool type, const FourVector &p) {
            Complex rpp(sqrt(PPlus(p))), rpm(sqrt(PMinus(p))), pt(PT(p));
            if( pt != Complex(0.0, 0.0))
                rpm = Complex(pt.real(), type ? pt.imag() : -pt.imag()) / rpp;

            m_type = type;
            m_u1 = rpp;
            m_u2 = rpm;
        }
        WeylSpinor(bool type, Complex u1, Complex u2) : m_type{type}, m_u1{u1}, m_u2{u2} {}
        WeylSpinor(const WeylSpinor&) = default;
        WeylSpinor(WeylSpinor&&) = default;
        WeylSpinor &operator=(const WeylSpinor&) = default;
        WeylSpinor &operator=(WeylSpinor&&) = default;
        ~WeylSpinor() = default;

        Complex U1() const { return m_u1; }
        Complex U2() const { return m_u2; }

        inline Complex operator*(const WeylSpinor &s) const {
            return m_u1*s.m_u2 - m_u2*s.m_u1;
        }

        WeylSpinor operator-() const {
            return {m_type, -m_u1, -m_u2};
        }

    private:
        static inline Complex PT(const FourVector &p) {
            return {p[1], p[2]};
        }

        static inline Complex PPlus(const FourVector &p) {
            return {p[0]+p[3], 0};
        }

        static inline Complex PMinus(const FourVector &p) {
            return {p[0]-p[3], 0};
        }

        bool m_type{};
        Complex m_u1, m_u2;
};

class SpinMatrix;

class Spinor {
    public:
        using Complex = std::complex<double>;

        Spinor() = default;
        Spinor(bool, bool, const int&, const FourVector &p, int ms=1);
        Spinor(bool type, bool bar, int hel, const std::array<Complex, 4> &u)
            : m_type{type}, m_bar{bar}, m_hel{hel}, m_u{u} {}
        Spinor(const Spinor&) = default;
        Spinor(Spinor&&) = default;
        Spinor& operator=(const Spinor&) = default;
        Spinor& operator=(Spinor&&) = default;

        const Complex& operator[](size_t i) const { return m_u[i]; }
        Complex& operator[](size_t i) { return m_u[i]; }
        Spinor Bar() const;
        bool Barred() const { return m_bar; }
        bool Type() const { return m_type; }
        int Helicity() const { return m_hel; }

        template<typename OStream>
        friend OStream& operator<<(OStream &os, const Spinor &spinor) {
            return os<< (!spinor.m_bar ?
                    (spinor.m_type ? "|u(" : "|v(") :
                    (spinor.m_type ? "<u(" : "<v("))
                   << (spinor.m_hel > 0 ? "+" : "-") << ") { "
                   << spinor[0] << ", " << spinor[1] << ", " << spinor[2] << " ," << spinor[3] << " }"
                   << (!spinor.m_bar ? ">" : "|");
        }

        Spinor operator-() const {
            return {m_type, m_bar, m_hel, {-m_u[0], -m_u[1], -m_u[2], -m_u[3]}};
        }

        Complex operator*(const Spinor &other) const {
            if(!m_bar) throw std::logic_error("Achilles::SpinorError: LHS spinor should be barred");
            if(other.m_bar) throw std::logic_error("Achilles::SpinorError: RHS spinor should not be barred");

            return m_u[0]*other[0] + m_u[1]*other[1] + m_u[2]*other[2] + m_u[3]*other[3];
        }

        SpinMatrix outer(const Spinor &other) const;

    private:
        bool m_type, m_bar;
        int m_hel;
        std::array<Complex, 4> m_u;
        static constexpr double tol = 1e-8;
};

inline Spinor UBarSpinor(int hel, const FourVector& mom) { return Spinor(true, true, hel, mom); }
inline Spinor VBarSpinor(int hel, const FourVector& mom) { return Spinor(false, true, hel, mom); }
inline Spinor USpinor(int hel, const FourVector& mom) { return Spinor(true, false, hel, mom); }
inline Spinor VSpinor(int hel, const FourVector& mom) { return Spinor(false, false, hel, mom); }

class SpinMatrix {
    public:
        using Complex = std::complex<double>;

        constexpr SpinMatrix() = default;
        constexpr SpinMatrix(std::array<Complex, 16> mat) : m_mat{std::move(mat)} {}

        static constexpr SpinMatrix Identity() {
            return {{1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1}};
        }

        static constexpr SpinMatrix Gamma_0() {
            return {{0, 0, 1, 0,
                     0, 0, 0, 1,
                     1, 0, 0, 0,
                     0, 1, 0, 0}};
        }

        static constexpr SpinMatrix Gamma_1() {
            return {{0, 0, 0, 1,
                     0, 0, 1, 0,
                     0, -1, 0, 0,
                     -1, 0, 0, 0}};
        }

        static constexpr SpinMatrix Gamma_2() {
            return {{0, 0, 0, Complex(0, -1),
                     0, 0, Complex(0, 1), 0,
                     0, Complex(0, 1), 0, 0,
                     Complex(0, -1), 0, 0, 0}};
        }

        static constexpr SpinMatrix Gamma_3() {
            return {{0, 0, 1, 0,
                     0, 0, 0, -1,
                     -1, 0, 0, 0,
                     0, 1, 0, 0}};
        }

        static constexpr SpinMatrix Gamma_5() {
            return {{-1, 0, 0, 0,
                     0, -1, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1}};
        }

        static SpinMatrix GammaMu(size_t i) {
            if(i == 0) return Gamma_0();
            else if(i == 1) return Gamma_1();
            else if(i == 2) return Gamma_2();
            else if(i == 3) return Gamma_3();
            else if(i == 5) return Gamma_5();
            else throw std::logic_error(fmt::format("Invalid Gamma Matrix: {}", i));
        }

        static SpinMatrix PL();
        static SpinMatrix PR();
        static SpinMatrix Slashed(const FourVector &mom);
        static SpinMatrix SigmaMuNu(size_t mu, size_t nu);

        SpinMatrix operator+=(const SpinMatrix &other) {
            for(size_t i = 0; i < m_mat.size(); ++i) {
                m_mat[i] += other[i];
            }
            return *this;
        }

        SpinMatrix operator+(const SpinMatrix &other) const {
            SpinMatrix result;
            for(size_t i = 0; i < result.size(); ++i) {
                result[i] = m_mat[i] + other[i];
            }

            return result;
        }

        SpinMatrix operator-(const SpinMatrix &other) const {
            SpinMatrix result;
            for(size_t i = 0; i < result.size(); ++i) {
                result[i] = m_mat[i] - other[i];
            }

            return result;
        }

        SpinMatrix operator-() const {
            SpinMatrix result;
            for(size_t i = 0; i < result.size(); ++i) {
                result[i] = -m_mat[i];
            }

            return result;
        }

        bool operator==(const SpinMatrix &other) const {
            return m_mat == other.m_mat;
        }

        Complex& operator[](size_t i) { return m_mat[i]; }
        const Complex& operator[](size_t i) const { return m_mat[i]; }

        constexpr size_t size() { return 16; }

        template<typename OStream>
        friend OStream& operator<<(OStream &os, const SpinMatrix &mat) {
            os << "SpinMatrix{";
            for(size_t i = 0; i < 4; ++i) {
                os << "{ ";
                for(size_t j = 0; j < 4; ++j) {
                    os << mat[4*i+j]; 
                    if(j != 3) os << ", ";
                }
                os << " }";
                if(i != 3) os << ", ";
            }
            os << "}";

            return os;
        }

    private:
        std::array<Complex, 16> m_mat{};
};

template<class T>
struct is_numeric 
    : std::integral_constant<
        bool,
        std::is_same<int, typename std::remove_cv<T>::type>::value ||
        std::is_same<float, typename std::remove_cv<T>::type>::value ||
        std::is_same<double, typename std::remove_cv<T>::type>::value ||
        std::is_same<long double, typename std::remove_cv<T>::type>::value ||
        std::is_same<std::complex<float>, typename std::remove_cv<T>::type>::value ||
        std::is_same<std::complex<double>, typename std::remove_cv<T>::type>::value ||
        std::is_same<std::complex<long double>, typename std::remove_cv<T>::type>::value> {};

template<typename T,
         std::enable_if_t<achilles::is_numeric<T>::value, bool> = true>
SpinMatrix operator*(const SpinMatrix &m, const T& scale) {
    SpinMatrix result;
    for(size_t i = 0; i < result.size(); ++i) {
        result[i] = m[i]*SpinMatrix::Complex(scale);
    }
    return result;
}

template<typename T,
         std::enable_if_t<achilles::is_numeric<T>::value, bool> = true>
SpinMatrix operator*(const T& scale, const SpinMatrix &m) {
    SpinMatrix result;
    for(size_t i = 0; i < result.size(); ++i) {
        result[i] = m[i]*SpinMatrix::Complex(scale);
    }
    return result;
}

inline SpinMatrix operator*(const SpinMatrix &lhs, const SpinMatrix &rhs) {
    SpinMatrix result;
    for(size_t i = 0; i < 4; ++i) {
        for(size_t j = 0; j < 4; ++j) {
            for(size_t k = 0; k < 4; ++k) {
                result[4*i + k] += lhs[4*i+j]*rhs[4*j+k];
            }
        }
    }

    return result;
}

Spinor operator*(const Spinor &lhs, const SpinMatrix &rhs);
Spinor operator*(const SpinMatrix &lhs, const Spinor &rhs);

template<typename T,
         std::enable_if_t<achilles::is_numeric<T>::value, bool> = true>
SpinMatrix operator/(const SpinMatrix &m, const T& scale) { 
    return m * (1.0/scale);
}

}

#endif
