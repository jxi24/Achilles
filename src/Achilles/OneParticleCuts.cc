#include "Achilles/OneParticleCuts.hh"
#include "Achilles/Constants.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Units.hh"

bool achilles::EnergyCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.E());
}

bool achilles::MomentumCut::MakeCut(const FourVector &mom) const {
    spdlog::trace("Momentum = {}, Cut Result = {}", mom.P(), CheckCut(mom.P()));
    return CheckCut(mom.P());
}

bool achilles::AngleThetaCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.Theta()/1.0_deg);
}

bool achilles::TransverseMomentumCut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.Pt());
}

bool achilles::ETheta2Cut::MakeCut(const FourVector &mom) const {
    return CheckCut(mom.E()*pow(mom.Theta(), 2));
}

bool achilles::Q0_QECut::MakeCut(const FourVector &mom) const {
    double q0 = (Constant::mN2-pow(Constant::mN-eb, 2)-pow(ParticleInfo(PID::muon()).Mass(), 2)+2*(mom.E()-mom.P()*mom.CosTheta())*mom.E())/(2*(Constant::mN-eb)-mom.E()+mom.P()*mom.CosTheta());
    return CheckCut(q0);
}
