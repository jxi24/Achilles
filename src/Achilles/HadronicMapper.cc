#include "Achilles/HadronicMapper.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/ParticleInfo.hh"
#include "spdlog/spdlog.h"

using achilles::QESpectralMapper;
using achilles::CoherentMapper;

void CoherentMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &) {
    point[HadronIdx()] = {ParticleInfo(PID::carbon()).Mass(), 0, 0, 0};
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, {});
}

double CoherentMapper::GenerateWeight(const std::vector<FourVector> &, std::vector<double>&) {
    return 1;
}

void QESpectralMapper::GeneratePoint(std::vector<FourVector> &point, const std::vector<double> &rans) {
    // Generate inital nucleon state
    double pmin = point[1].E() - sqrt(pow(point[1].E(), 2) + 2*point[1].E()*Constant::mN + Constant::mN2 - Smin());
    double pmax = point[1].E() + sqrt(pow(point[1].E(), 2) + 2*point[1].E()*Constant::mN + Constant::mN2 - Smin());
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    const double mom = dp*rans[0] + pmin;
    double cosT_max = (2*point[1].E()*Constant::mN+Constant::mN2-mom*mom-Smin())/(2*point[1].E()*mom);
    cosT_max = cosT_max > 1 ? 1 : cosT_max;
    const double cosT = (cosT_max + 1)*rans[1] - 1;
    const double sinT = sqrt(1 - cosT*cosT);
    const double phi = dPhi*rans[2];
    ThreeVector pmom = {mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};

    const double det = pow(point[1].E(), 2) + mom*mom + 2*pmom*point[1].Vec3() + Smin();
    double emax = Constant::mN + point[1].E() - sqrt(det);
    emax = emax > 400 ? 400 : emax;
    const double energy = emax*rans[3] - 1e-8;
    // if(emax < 0) energy = emax - 1;
    // const double energy = dE*rans[3];
    
    // double cosT_max = (Constant::mN2 + energy*energy - 2*Constant::mN*energy - mom*mom + 2*point[1].E()*Constant::mN - 2*point[1].E()*energy - Smin())/(2*mom*point[1].P());
    // cosT_max = cosT_max > 1 ? 1 : cosT_max;
    // const double cosT = (cosT_max + 1) * rans[1] - 1;
    // const double sinT = sqrt(1 - cosT*cosT);

    point[HadronIdx()] = {Constant::mN - energy, mom*sinT*cos(phi), mom*sinT*sin(phi), mom*cosT};
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  dp = {}", dp);
    spdlog::trace("  cosT_min = {}", cosT_max);
    spdlog::trace("  cosT_max = {}", cosT_max);
    spdlog::trace("  cosT = {}", cosT);
    spdlog::trace("  mom = {}", mom);
    spdlog::trace("  energy = {}", energy);
    spdlog::trace("  s = {}", (point[0] + point[1]).M2());
    spdlog::trace("  s_min = {}", Smin());
}

double QESpectralMapper::GenerateWeight(const std::vector<FourVector> &point, std::vector<double> &rans) {
    double pmin = point[1].E() - sqrt(pow(point[1].E(), 2) + 2*point[1].E()*Constant::mN + Constant::mN2 - Smin());
    double pmax = point[1].E() + sqrt(pow(point[1].E(), 2) + 2*point[1].E()*Constant::mN + Constant::mN2 - Smin());
    pmin = pmin < 0 ? 0 : pmin;
    pmax = pmax > 800 ? 800 : pmax;
    double dp = pmax - pmin;
    rans[0] = (point[HadronIdx()].P() - pmin)/dp;
    double cosT_max = (2*point[1].E()*Constant::mN+Constant::mN2-point[0].P2()-Smin())/(2*point[1].E()*point[0].P());
    cosT_max = cosT_max > 1 ? 1 : cosT_max;
    double dCos = (cosT_max + 1);
    rans[1] = (point[HadronIdx()].CosTheta()+1)/dCos;
    rans[2] = point[HadronIdx()].Phi()/dPhi;

    const double det = pow(point[1].E(), 2) + point[0].P2() 
                     + 2*point[0].Vec3()*point[1].Vec3() + Smin();
    double emax = Constant::mN + point[1].E() - sqrt(det);
    emax = emax > 400 ? 400 : emax;
    const double energy = Constant::mN - point[HadronIdx()].E();
    // if(energy < 0) return std::numeric_limits<double>::infinity();
    const double dE = emax;
    rans[3] = (energy + 1e-8)/emax;
    // rans[3] = (Constant::mN - point[HadronIdx()].E())/dE; 

    // double cosT_max = (point[HadronIdx()].M2() + 2*point[1].E()*point[HadronIdx()].E() - Smin())/(2*point[HadronIdx()].P()*point[1].P());
    // cosT_max = cosT_max > 1 ? 1 : cosT_max;
    // const double dCos = (cosT_max + 1);
    rans[1] = (point[HadronIdx()].CosTheta()+1)/dCos;

    double wgt = 1.0/point[0].P2()/dp/dCos/dPhi/dE;
    Mapper<FourVector>::Print(__PRETTY_FUNCTION__, point, rans);
    spdlog::trace("  Weight: {}", wgt);
    spdlog::trace("  dp: {}", dp);
    spdlog::trace("  dCos: {}", dCos);
    spdlog::trace("  dPhi: {}", dPhi);
    spdlog::trace("  dE: {}", dE);

    return wgt;
}
