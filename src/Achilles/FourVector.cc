#include <iostream>

#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"
#include "spdlog/spdlog.h"

using achilles::FourVector;
using achilles::ThreeVector;

FourVector::FourVector(const ThreeVector& other, const double& E) noexcept
    : vec({{E, other[0], other[1], other[2]}}) {}

double FourVector::M() const noexcept {
    return std::sqrt(M2());
}

double FourVector::Theta() const noexcept {
    return atan2(Pt(), Pz());
}

double FourVector::Phi() const noexcept {
    const double phi = atan2(Py(), Px());
    if(phi < 0) return phi + 2*M_PI;
    return phi;
}

double FourVector::Rapidity() const noexcept {
    return log((E() + Pz())/(E() - Pz()))/2;
}

double FourVector::DeltaR(const FourVector& other) const noexcept {
    const double DEta = Rapidity() - other.Rapidity();
    const double DPhi = Phi() - other.Phi();
    return sqrt(DEta*DEta + DPhi*DPhi);
}

double FourVector::CosAngle(const FourVector &other) const noexcept {
    auto p1 = Vec3(), p2 = other.Vec3();
    return std::max(std::min(p1*p2/(p1.P()*p2.P()), 1.0), -1.0);
}

double FourVector::Angle(const FourVector &other) const noexcept {
    return std::acos(CosAngle(other));
}

ThreeVector FourVector::Vec3() const noexcept {
    return {Px(), Py(), Pz()};
}

void FourVector::SetVectM(const ThreeVector& vec3, const double& mass) noexcept {
    Px() = vec3.Px();
    Py() = vec3.Py();
    Pz() = vec3.Pz();
    E() = sqrt(mass*mass + vec3*vec3);
}

FourVector FourVector::Boost(const ThreeVector& beta) const noexcept {
    const double beta2 = beta*beta;
    const double gamma = 1.0/sqrt(1.0 - beta2);
    const double betap = beta[0]*Px() + beta[1]*Py() + beta[2]*Pz();
    const double gamma2 = beta2 > 0 ? (gamma-1.0)/beta2 : 0.0;

    double pX = Px() + gamma2*betap*beta[0]+gamma*beta[0]*E();
    double pY = Py() + gamma2*betap*beta[1]+gamma*beta[1]*E();
    double pZ = Pz() + gamma2*betap*beta[2]+gamma*beta[2]*E();
    double energy = gamma*(E() + betap);

    return {energy, pX, pY, pZ};
}

FourVector FourVector::Boost(const double& beta_x, const double& beta_y,
                                             const double& beta_z) const noexcept {
    return Boost(ThreeVector(beta_x, beta_y, beta_z));
}

double FourVector::SmallOMCT(const FourVector& v) const noexcept {
  double mag(sqrt(P2()*v.P2()));
  double pq(vec[1]*v[1]+vec[2]*v[2]+vec[3]*v[3]);
  double ct(std::min(std::max(pq/mag,-1.),1.));
  if (ct<0.) return 1.-ct;
  double st(this->Vec3().Cross(v.Vec3()).P()/mag);
  double st2(st/(2.*sqrt((ct+1)/2.)));
  return 2.*st2*st2;
}

double FourVector::SmallMLDP(const FourVector& v) const noexcept {
  return vec[0]*v[0]*SmallOMCT(v);
}

FourVector FourVector::Rotate(const RotMat &mat) const noexcept {
    return {E(),
            mat[0]*Px()+mat[1]*Py()+mat[2]*Pz(),
            mat[3]*Px()+mat[4]*Py()+mat[5]*Pz(),
            mat[6]*Px()+mat[7]*Py()+mat[8]*Pz()};
}

FourVector FourVector::RotateBack(const RotMat &mat) const noexcept {
    return {E(),
            mat[0]*Px()+mat[3]*Py()+mat[6]*Pz(),
            mat[1]*Px()+mat[4]*Py()+mat[7]*Pz(),
            mat[2]*Px()+mat[5]*Py()+mat[8]*Pz()};
}

achilles::FourVector::RotMat FourVector::Align(const ThreeVector &axis) const noexcept {

    ThreeVector a = Vec3().Unit();

    auto v = a.Cross(axis);
    double c = a.Dot(axis);

    return {{1-v[1]*v[1]/(1+c)-v[2]*v[2]/(1+c), -v[2]+v[0]*v[1]/(1+c), v[1]+v[0]*v[2]/(1+c),
            v[2]+v[0]*v[1]/(1+c), 1-v[0]*v[0]/(1+c)-v[2]*v[2]/(1+c), -v[0]+v[1]*v[2]/(1+c),
            -v[1]+v[0]*v[2]/(1+c), v[0]+v[1]*v[2]/(1+c), 1-v[0]*v[0]/(1+c)-v[1]*v[1]/(1+c)}};
}

achilles::FourVector::RotMat FourVector::AlignZ() const noexcept {
    ThreeVector z{0, 0, 1};
    return Align(z);
}

FourVector FourVector::Cross(const FourVector& other) const noexcept {
    ThreeVector result = this -> Vec3().Cross(other.Vec3());
    return {result, 0};
}

ThreeVector FourVector::BoostVector() const noexcept {
    return {Px()/E(), Py()/E(), Pz()/E()};
}

FourVector& FourVector::operator+=(const FourVector& other) noexcept {
    E() += other.E();
    Px() += other.Px();
    Py() += other.Py();
    Pz() += other.Pz();

    return *this;
}

FourVector& FourVector::operator-=(const FourVector& other) noexcept {
    E() -= other.E();
    Px() -= other.Px();
    Py() -= other.Py();
    Pz() -= other.Pz();

    return *this;
}

FourVector& FourVector::operator*=(const double& scale) noexcept {
    E() *= scale;
    Px() *= scale;
    Py() *= scale;
    Pz() *= scale;

    return *this;
}

FourVector& FourVector::operator/=(const double& scale) {
    E() /= scale;
    Px() /= scale;
    Py() /= scale;
    Pz() /= scale;

    return *this;
}

double FourVector::operator*(const FourVector& other) const noexcept {
    return E()*other.E()
        - (Px()*other.Px() + Py()*other.Py() + Pz()*other.Pz());
}

FourVector FourVector::operator-() const noexcept {
   return {-vec[0],-vec[1],-vec[2],-vec[3]};
}

FourVector FourVector::operator+() const noexcept {
    return *this;
}

FourVector FourVector::operator*(const double& scale) const noexcept {
    return FourVector(*this)*=scale;
}

FourVector FourVector::operator/(const double& scale) const {
    return FourVector(*this)/=scale;
}

FourVector FourVector::operator+(const FourVector& other) const noexcept {
    return FourVector(*this)+=other;
}

FourVector FourVector::operator-(const FourVector& other) const noexcept {
    return FourVector(*this)-=other;
}

bool FourVector::operator==(const FourVector& other) const noexcept {
    return vec==other.vec;
}

std::string FourVector::ToString() const noexcept {
    return "FourVector(" + std::to_string(vec[0]) + ", " + std::to_string(vec[1])
        + ", " + std::to_string(vec[2]) + ", " + std::to_string(vec[3]) + ")";
}

namespace achilles {

std::istream& operator>>(std::istream& is, FourVector& vec) {
    std::string head_name = "FourVector(";
    std::string head(head_name.size(), ' '), sep1(1, ' '), sep2(1, ' '),
        sep3(1, ' '), tail(1, ' ');
    double e{}, px{}, py{}, pz{};
    is.read(&head[0], 11);
    is >> e;
    is.read(&sep1[0], 1);
    is >> px;
    is.read(&sep2[0], 1);
    is >> py;
    is.read(&sep3[0], 1);
    is >> pz;
    is.read(&tail[0], 1);
    if(head == head_name &&
       sep1 == "," && sep2 == "," && 
       sep3 == "," && tail == ")") 
        vec = FourVector(e, px, py, pz);
    return is;
}

FourVector operator*(const double& s, const FourVector& v) noexcept {
    return v*s;
}

}
