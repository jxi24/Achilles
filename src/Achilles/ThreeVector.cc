#include <iostream>

#include "Achilles/ThreeVector.hh"

using achilles::ThreeVector;

double ThreeVector::Theta() const noexcept {
    return atan2(Pt(), Pz());
}

double ThreeVector::Phi() const noexcept {
    const double phi = atan2(Py(), Px());
    if(phi < 0) return phi + 2*M_PI;
    return phi;
}

ThreeVector ThreeVector::Cross(const ThreeVector& other) const noexcept {
    const double px = vec[1]*other.vec[2] - vec[2]*other.vec[1];    
    const double py = vec[2]*other.vec[0] - vec[0]*other.vec[2];
    const double pz = vec[0]*other.vec[1] - vec[1]*other.vec[0];    

    return {px, py, pz};
}

ThreeVector ThreeVector::Unit() const noexcept {
    const double norm = Magnitude();

    return {{vec[0]/norm, vec[1]/norm, vec[2]/norm}};
}

ThreeVector ThreeVector::Rotate(const std::array<double, 3> &angles) const noexcept {
    const double c1 = cos(angles[0]), s1 = sin(angles[0]);
    const double c2 = cos(angles[1]), s2 = sin(angles[1]);
    const double c3 = cos(angles[2]), s3 = sin(angles[2]);

    return {{(c1*c3-c2*s1*s3)*vec[0]+(-c1*s3-c2*c3*s1)*vec[1]+s1*s2*vec[2],
            (c3*s1+c1*c2*s3)*vec[0]+(c1*c2*c3-s1*s3)*vec[1]-c1*s2*vec[2],
            s2*s3*vec[0]+c3*s2*vec[1]+c2*vec[2]}}; 
}

ThreeVector ThreeVector::Rotate(const RotMat &mat) const noexcept {
    return {{mat[0]*vec[0]+mat[1]*vec[1]+mat[2]*vec[2],
            mat[3]*vec[0]+mat[4]*vec[1]+mat[5]*vec[2],
            mat[6]*vec[0]+mat[7]*vec[1]+mat[8]*vec[2]}};
}

achilles::ThreeVector::RotMat ThreeVector::Align(const ThreeVector &axis) const noexcept {

    ThreeVector a = Unit();

    auto v = a.Cross(axis);
    double c = a.Dot(axis);

    return {{1-v[1]*v[1]/(1+c)-v[2]*v[2]/(1+c), -v[2]+v[0]*v[1]/(1+c), v[1]+v[0]*v[2]/(1+c),
            v[2]+v[0]*v[1]/(1+c), 1-v[0]*v[0]/(1+c)-v[2]*v[2]/(1+c), -v[0]+v[1]*v[2]/(1+c),
            -v[1]+v[0]*v[2]/(1+c), v[0]+v[1]*v[2]/(1+c), 1-v[0]*v[0]/(1+c)-v[1]*v[1]/(1+c)}};
}

achilles::ThreeVector::RotMat ThreeVector::AlignZ() const noexcept {
    ThreeVector z{0, 0, 1};
    return Align(z);
}

ThreeVector& ThreeVector::operator+=(const ThreeVector& other) noexcept {
    vec[0] += other.vec[0];
    vec[1] += other.vec[1];
    vec[2] += other.vec[2];

    return *this;
}

ThreeVector& ThreeVector::operator-=(const ThreeVector& other) noexcept {
    vec[0] -= other.vec[0];
    vec[1] -= other.vec[1];
    vec[2] -= other.vec[2];

    return *this;
}

ThreeVector& ThreeVector::operator*=(const double& scale) noexcept {
    vec[0] *= scale;
    vec[1] *= scale;
    vec[2] *= scale;

    return *this;
}

ThreeVector& ThreeVector::operator/=(const double& scale) {
    vec[0] /= scale;
    vec[1] /= scale;
    vec[2] /= scale;

    return *this;
}

double ThreeVector::operator*(const ThreeVector& other) const noexcept {
    return vec[0]*other.vec[0] + vec[1]*other.vec[1] + vec[2]*other.vec[2];
}

ThreeVector ThreeVector::operator-() const noexcept {
   return {-vec[0], -vec[1], -vec[2]};
}

ThreeVector ThreeVector::operator+() const noexcept {
    return *this;
}

ThreeVector ThreeVector::operator*(const double& scale) const noexcept {
    return ThreeVector(*this)*=scale;
}

ThreeVector ThreeVector::operator/(const double& scale) const {
    return ThreeVector(*this)/=scale;
}

ThreeVector ThreeVector::operator+(const ThreeVector& other) const noexcept {
    return ThreeVector(*this)+=other;
}

ThreeVector ThreeVector::operator-(const ThreeVector& other) const noexcept {
    return ThreeVector(*this)-=other;
}

bool ThreeVector::operator==(const ThreeVector& other) const noexcept {
    return vec==other.vec;
}

std::string ThreeVector::ToString() const noexcept {
    return "ThreeVector(" + std::to_string(vec[0]) + ", " 
        + std::to_string(vec[1]) + ", " + std::to_string(vec[2]) + ")";
}

namespace achilles {

std::istream& operator>>(std::istream& is, ThreeVector& vec) {
    std::string head(12, ' '), sep1(1, ' '), sep2(1, ' '), tail(1, ' ');
    double px, py, pz;
    is.read(&head[0], 12);
    is >> px;
    is.read(&sep1[0], 1);
    is >> py;
    is.read(&sep2[0], 1);
    is >> pz;
    is.read(&tail[0], 1);
    if(head == "ThreeVector(" &&
       sep1 == "," && sep2 == "," && 
       tail == ")") 
        vec = ThreeVector(px, py, pz);
    return is;
}

ThreeVector operator*(const double& s, const ThreeVector& v) noexcept {
    return v*s;
}

}
