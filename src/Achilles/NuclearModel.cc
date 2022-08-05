#include "Achilles/NuclearModel.hh"
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Spinor.hh"
#include "Achilles/Particle.hh"

using achilles::NuclearModel;
using achilles::Coherent;
using achilles::QESpectral;
using Type = achilles::FormFactorInfo::Type;

NuclearModel::NuclearModel(const YAML::Node& config,
                           FormFactorBuilder &ffbuilder = FormFactorBuilder::Instance()) {
    spdlog::debug("Setting up form factors");
    const auto vectorFF = config["vector"].as<std::string>();
    const auto axialFF = config["axial"].as<std::string>();
    const auto coherentFF = config["coherent"].as<std::string>();
    m_form_factor = ffbuilder.Vector(vectorFF, config[vectorFF])
                             .AxialVector(axialFF, config[axialFF])
                             .Coherent(coherentFF, config[coherentFF])
                             .build();
}

NuclearModel::FormFactorMap NuclearModel::CouplingsFF(const FormFactor::Values &formFactors,
                                                      const std::vector<FormFactorInfo> &ffInfo) const {
    FormFactorMap results{};

    for(const auto & ff : ffInfo) {
        spdlog::trace("Form Factor: {}, Coupling: {}", ff.form_factor, ff.coupling);
        switch(ff.form_factor) {
            case Type::F1p:
                results[Type::F1] += formFactors.F1p*ff.coupling;
                break;
            case Type::F1n:
                results[Type::F1] += formFactors.F1n*ff.coupling;
                break;
            case Type::F2p:
                results[Type::F2] += formFactors.F2p*ff.coupling;
                break;
            case Type::F2n:
                results[Type::F2] += formFactors.F2n*ff.coupling;
                break;
            case Type::FA:
                results[Type::FA] += formFactors.FA*ff.coupling;
                break;
            case Type::FCoh:
                results[Type::FCoh] += formFactors.Fcoh*ff.coupling;
                break;
            case Type::F1:
            case Type::F2:
                throw std::runtime_error("Types F1 and F2 are reserved for nuclear models");
        }
    }

    return results;
}

YAML::Node NuclearModel::LoadFormFactor(const YAML::Node &config) {
    return YAML::LoadFile(config["NuclearModel"]["FormFactorFile"].as<std::string>());
}

// TODO: Clean this up such that the nucleus isn't loaded twice
Coherent::Coherent(const YAML::Node &config, const YAML::Node &form_factor,
                   FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, builder) {
    nucleus_pid = config["Nucleus"].as<Nucleus>().ID();
}

std::vector<achilles::NuclearModel::Currents> Coherent::CalcCurrents(const Event &event,
                                                                   const std::vector<FFInfoMap> &ff) const {
    auto pIn = event.Momentum().front();
    auto pOut = event.Momentum()[2];
    auto qVec = event.Momentum()[1];
    for(size_t  i = 3; i < event.Momentum().size(); ++i) {
        qVec -= event.Momentum()[i];
    }

    auto ffVals = EvalFormFactor(-qVec.M2());

    // Calculate coherent contributions
    std::vector<Currents> results(1);
    for(const auto &formFactor : ff[2]) {
        auto ffVal = CouplingsFF(ffVals, formFactor.second);
        spdlog::trace("fcoh = {}", ffVal[Type::FCoh]);

        Current current;
        std::vector<std::complex<double>> subcur(4);
        for(size_t i = 0; i < subcur.size(); ++i) {
            subcur[i] = (pIn[i] + pOut[i])*ffVal[Type::FCoh];
        }
        current.push_back(subcur);
        results[0][formFactor.first] = current;
        spdlog::trace("HadronicCurrent[{}] = [{}, {}, {}, {}]", formFactor.first,
                      subcur[0], subcur[1], subcur[2], subcur[3]);
    }

    return results;
}

void Coherent::AllowedStates(Process_Info &info) const {
    // Check for charge conservation
    int charge = -ParticleInfo(info.m_ids[0]).IntCharge();
    for(size_t i = 1; i < info.m_ids.size(); ++i) {
        charge += ParticleInfo(info.m_ids[i]).IntCharge();
    }
    charge /= 3;
    spdlog::debug("Charge = {}", charge);
    if(charge != 0)
        throw std::runtime_error(fmt::format("Coherent: Requires charge 0, but found charge {}", charge));

    info.m_states[{nucleus_pid}] = {nucleus_pid}; 
}

bool Coherent::FillNucleus(Event &event, const std::vector<double> &xsecs) const {
    // Calculate total cross section
    if(xsecs[0] == 0) return false;
    event.SetMEWeight(xsecs[0]);

    // Remove all nucleons
    event.CurrentNucleus() -> Nucleons().clear();

    // Setup initial and final state nucleus
    Particle initial = Particle(nucleus_pid, event.Momentum().front());
    initial.Status() = ParticleStatus::initial_state;
    event.CurrentNucleus() -> Nucleons().push_back(initial);
    Particle final(nucleus_pid, event.Momentum()[2]);
    final.Status() = ParticleStatus::final_state;
    event.CurrentNucleus() -> Nucleons().push_back(final);

    return true;
}

std::unique_ptr<NuclearModel> Coherent::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<Coherent>(config, form_factor);
}

// TODO: Clean this interface up
QESpectral::QESpectral(const YAML::Node &config, const YAML::Node &form_factor,
                       FormFactorBuilder &builder = FormFactorBuilder::Instance()) 
        : NuclearModel(form_factor, builder),
          spectral_proton{config["NuclearModel"]["SpectralP"].as<std::string>()},
          spectral_neutron{config["NuclearModel"]["SpectralN"].as<std::string>()} {
    b_ward = config["NuclearModel"]["Ward"].as<bool>();
}

std::vector<NuclearModel::Currents> QESpectral::CalcCurrents(const Event &event,
                                                             const std::vector<FFInfoMap> &ff) const {
    auto pIn = event.Momentum().front();
    auto pOut = event.Momentum()[2];
    auto qVec = event.Momentum()[1];
    for(size_t i = 3; i < event.Momentum().size(); ++i) {
        qVec -= event.Momentum()[i];
    }
    auto removal_energy = Constant::mN - pIn.E();
    auto free_energy = sqrt(pIn.P2() + Constant::mN2);
    auto ffVals = EvalFormFactor(-qVec.M2()/1.0_GeV/1.0_GeV);
    auto omega = qVec.E();
    qVec.E() = qVec.E() + pIn.E() - free_energy;

    std::vector<Currents> results(2);
    // TODO: Move to the phase space definition
    std::vector<double> spectral(2);
    spectral[0] = spectral_proton(pIn.P(), removal_energy);
    spectral[1] = spectral_neutron(pIn.P(), removal_energy);
    spdlog::debug("Spectral function: S_p({}, {}) = {}, S_n({}, {}) = {}",
                  removal_energy, pIn.P(), spectral[0],
                  removal_energy, pIn.P(), spectral[1]);
    // Setup spinors
    pIn.E() = free_energy;
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pOut);
    ubar[1] = UBarSpinor(1, pOut);
    u[0] = USpinor(-1, -pIn);
    u[1] = USpinor(1, -pIn);

    // Loop over proton and neutron
    for(size_t i = 0; i < results.size(); ++i) {
        // Calculate nucleon contributions
        for(const auto &formFactor : ff[i]) {
            auto ffVal = CouplingsFF(ffVals, formFactor.second);
            spdlog::debug("{}: f1 = {}, f2 = {}, fa = {}",
                          i, ffVal[Type::F1], ffVal[Type::F2], ffVal[Type::FA]);
            auto current = HadronicCurrent(ubar, u, qVec, ffVal);
            for(auto &subcur : current) {
                for(auto &val : subcur) {
                    // TODO: Move this to phase space 
                    // TODO: Move normalization to spectral function definition
                    val *= sqrt(spectral[i]/6);
                }
                // Correct the Ward identity
                if(b_ward) subcur[3] = omega/qVec.P()*subcur[0];
            }
            results[i][formFactor.first] = current;
        }
    }
    return results;
}

void QESpectral::AllowedStates(Process_Info &info) const {
    // Check for charge conservation
    int charge = -ParticleInfo(info.m_ids[0]).IntCharge();
    for(size_t i = 1; i < info.m_ids.size(); ++i) {
        charge += ParticleInfo(info.m_ids[i]).IntCharge();
    }
    charge /= 3;
    if(std::abs(charge) > 1)
        throw std::runtime_error(fmt::format("Quasielastic: Requires |charge| < 2, but found |charge| {}", std::abs(charge)));

    switch(charge) {
        case -1: // Final state has less charge than initial
            info.m_states[{PID::neutron()}] = {PID::proton()}; 
            break; 
        case 0: // Same charge in inital and final
            info.m_states[{PID::proton()}] = {PID::proton()}; 
            info.m_states[{PID::neutron()}] = {PID::neutron()}; 
            break;
        case 1: // Final state has more charge than initial
            info.m_states[{PID::proton()}] = {PID::neutron()}; 
            break;
    }
}

bool QESpectral::FillNucleus(Event &event, const std::vector<double> &xsecs) const {
    // Calculate total cross section
    for(size_t i = 0; i < event.CurrentNucleus() -> Nucleons().size(); ++i) {
        auto current_nucleon = event.CurrentNucleus() -> Nucleons()[i];
        if(current_nucleon.ID() == PID::proton()) {
            event.MatrixElementWgt(i) = xsecs[0];
        } else {
            event.MatrixElementWgt(i) = xsecs[1];
        }
    }
    if(!event.TotalCrossSection())
        return false;

    return true;
}

std::unique_ptr<NuclearModel> QESpectral::Construct(const YAML::Node &config) {
    auto form_factor = LoadFormFactor(config);
    return std::make_unique<QESpectral>(config, form_factor);
}

NuclearModel::Current QESpectral::HadronicCurrent(const std::array<Spinor, 2> &ubar,
                                                  const std::array<Spinor, 2> &u,
                                                  const FourVector &qVec,
                                                  const FormFactorMap &ffVal) const {
    Current result;
    std::array<SpinMatrix, 4> gamma{};
    auto mpi2 = pow(ParticleInfo(211).Mass(), 2);
    auto ffAP = 2.0*Constant::mN2/(-qVec.M2()+mpi2)*ffVal.at(Type::FA);
    for(size_t mu = 0; mu < 4; ++mu) {
        gamma[mu] = ffVal.at(Type::F1)*SpinMatrix::GammaMu(mu);
        gamma[mu] += ffVal.at(Type::FA)*SpinMatrix::GammaMu(mu)*SpinMatrix::Gamma_5();
        // REMOVE: 0.0
        gamma[mu] += ffAP*SpinMatrix::Gamma_5()*qVec[mu]/Constant::mN;
        double sign = 1;
        for(size_t nu = 0; nu < 4; ++nu) {
            gamma[mu] += std::complex<double>(0, 1)*(ffVal.at(Type::F2)*SpinMatrix::SigmaMuNu(mu, nu)*sign*qVec[nu]/(2*Constant::mN));
            sign = -1;
        }
    }

    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            std::vector<std::complex<double>> subcur(4);
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i]*gamma[mu]*u[j];
            }
            result.push_back(subcur);
        }
    }

    return result;
}
