#include <iostream>
#include <utility>
#include  <fstream>
#include <cmath>

#include "Achilles/HardScattering.hh"
#include "Achilles/Constants.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/Event.hh"
#include "Achilles/Random.hh"

// Aliases for most common types
using achilles::Particles;
using achilles::HardScattering;
using achilles::LeptonicCurrent;

void LeptonicCurrent::Initialize(const Process_Info &process) {
    using namespace achilles::Constant;
    const std::complex<double> i(0, 1);
    // Determine process
    bool init_neutrino = ParticleInfo(process.m_ids[0]).IsNeutrino();
    bool neutral_current = NeutralCurrent(process.m_ids[0], process.m_ids[1]);
    bool charged_current = ChargedCurrent(init_neutrino, process.m_ids[0], process.m_ids[1]);
    if(!neutral_current && !charged_current) 
        throw std::runtime_error("HardScattering: Invalid process");

    // TODO: Define couplings correctly
    if(charged_current) {
        pid = init_neutrino ? (process.m_ids[0].AsInt() < 0 ? -24 : 24)
                            : (process.m_ids[0].AsInt() < 0 ? 24 : -24);
        coupl_right = 0;
        coupl_left = ee*i/(sw*sqrt(2));
        mass = Constant::MW;
        width = Constant::GAMW;
    } else if(neutral_current) {
        if(init_neutrino) {
            coupl_left = (cw*ee*i)/(2*sw)+(ee*i*sw)/(2*cw);
            coupl_right = 0;
            pid = 23;
            mass = Constant::MZ;
            width = Constant::GAMZ;
        } else {
            coupl_right = -ee*i;
            coupl_left = coupl_right;
            pid = 22;
        }
    }
    anti = process.m_ids[0].AsInt() < 0;
}

bool LeptonicCurrent::NeutralCurrent(achilles::PID initial, achilles::PID final) const {
    return initial == final;
}

bool LeptonicCurrent::ChargedCurrent(bool neutrino, achilles::PID initial, achilles::PID final) const {
    return initial.AsInt() - (2*neutrino - 1) == final.AsInt();
}

achilles::FFDictionary LeptonicCurrent::GetFormFactor() {
    FFDictionary results;
    static constexpr std::complex<double> i(0, 1);
    using namespace achilles::Constant;
    // TODO: Double check form factors
    if(pid == 24) {
        const std::complex<double> coupl = ee*i/(sw*sqrt(2)*2);
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F1n, -coupl},
                                         {FormFactorInfo::Type::F2p, coupl},
                                         {FormFactorInfo::Type::F2n, -coupl},
                                         {FormFactorInfo::Type::FA, coupl}};
        results[{PID::neutron(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == -24) {
        const std::complex<double> coupl = ee*i/(sw*sqrt(2)*2);
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                          {FormFactorInfo::Type::F1n, -coupl},
                                          {FormFactorInfo::Type::F2p, coupl},
                                          {FormFactorInfo::Type::F2n, -coupl},
                                          {FormFactorInfo::Type::FA, coupl}};
        results[{PID::proton(), pid}] = {};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 23) {
        const std::complex<double> coupl1 = cw*ee*i/(2*sw)-ee*i*sw/(2*cw);
        const std::complex<double> coupl2 = -(cw*ee*i/(2*sw));
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl1},
                                         {FormFactorInfo::Type::F1n, coupl2},
                                         {FormFactorInfo::Type::F2p, coupl1},
                                         {FormFactorInfo::Type::F2n, coupl2},
                                         {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl1},
                                          {FormFactorInfo::Type::F1p, coupl2},
                                          {FormFactorInfo::Type::F2n, coupl1},
                                          {FormFactorInfo::Type::F2p, coupl2},
                                          {FormFactorInfo::Type::FA, coupl2}};
        results[{PID::carbon(), pid}] = {};
    } else if(pid == 22) {
        const std::complex<double> coupl = i*ee;
        results[{PID::proton(), pid}] = {{FormFactorInfo::Type::F1p, coupl},
                                         {FormFactorInfo::Type::F2p, coupl}};
        results[{PID::neutron(), pid}] = {{FormFactorInfo::Type::F1n, coupl},
                                          {FormFactorInfo::Type::F2n, coupl}};
        results[{PID::carbon(), pid}] = {{FormFactorInfo::Type::FCoh, 6.0*coupl}};
    } else {
        throw std::runtime_error("LeptonicCurrent: Invalid probe");
    }

    return results;
}

achilles::Currents LeptonicCurrent::CalcCurrents(const std::vector<FourVector> &p,
                                                 const double&) const {
    Currents currents;

    // Setup spinors
    FourVector pU, pUBar;
    if(anti) {
        pUBar = -p[1];
        pU = p.back();
    } else {
        pU = -p[1];
        pUBar = p.back();
    }
    std::array<Spinor, 2> ubar, u;
    ubar[0] = UBarSpinor(-1, pUBar);
    ubar[1] = UBarSpinor(1, pUBar);
    u[0] = USpinor(-1, pU);
    u[1] = USpinor(1, pU);

    // Calculate currents
    Current result;
    double q2 = (p[1] - p.back()).M2();
    std::complex<double> prop = std::complex<double>(0, 1)/(q2-mass*mass-std::complex<double>(0, 1)*mass*width);
    spdlog::trace("Calculating Current for {}", pid);
    for(size_t i = 0; i < 2; ++i) {
        for(size_t j = 0; j < 2; ++j) {
            std::vector<std::complex<double>> subcur(4);
            for(size_t mu = 0; mu < 4; ++mu) {
                subcur[mu] = ubar[i]*(coupl_left*SpinMatrix::GammaMu(mu)*SpinMatrix::PL()
                                    + coupl_right*SpinMatrix::GammaMu(mu)*SpinMatrix::PR())*u[j]*prop;
                spdlog::trace("Current[{}][{}] = {}", 2*i+j, mu, subcur[mu]);
            }
            result.push_back(subcur);
        }
    }
    currents[pid] = result;

    return currents;
}

void HardScattering::SetProcess(const Process_Info &process) {
    spdlog::debug("Adding Process: {}", process);
    m_leptonicProcess = process;
    #ifndef ENABLE_BSM
    m_current.Initialize(process);
    SMFormFactor = m_current.GetFormFactor();
    #endif
}

achilles::Currents HardScattering::LeptonicCurrents(const std::vector<FourVector> &p,
                                                    const double &mu2) const {
#ifdef ENABLE_BSM
    // TODO: Move adapter code into Sherpa interface code
    std::vector<std::array<double, 4>> mom(p.size());
    std::vector<int> pids;
    for(const auto &elm : m_leptonicProcess.m_mom_map) {
        pids.push_back(static_cast<int>(elm.second));
        mom[elm.first] = (p[elm.first]/1_GeV).Momentum();
        spdlog::debug("PID: {}, Momentum: ({}, {}, {}, {})", pids.back(),
                      mom[elm.first][0], mom[elm.first][1], mom[elm.first][2], mom[elm.first][3]); 
    }
    auto currents = p_sherpa -> Calc(pids, mom, mu2);

    for(auto &current : currents) { 
        spdlog::trace("Current for {}", current.first);
        for(size_t i = 0; i < current.second.size(); ++i) {
            for(size_t j = 0; j < current.second[0].size(); ++j) {
                current.second[i][j] /= pow(1_GeV, static_cast<double>(mom.size())-3);
                spdlog::trace("Current[{}][{}] = {}", i, j, current.second[i][j]);
            }
        }
    }
    return currents;
#else
    return m_current.CalcCurrents(p, mu2);
#endif
}

std::vector<double> HardScattering::CrossSection(Event &event) const {
    // Calculate leptonic currents
    auto leptonCurrent = LeptonicCurrents(event.Momentum(), 100);

    // Calculate the hadronic currents
    // TODO: Clean this up and make generic for the nuclear model
    // TODO: Move this to initialization to remove check each time
    static std::vector<NuclearModel::FFInfoMap> ffInfo;
    if(ffInfo.empty()) {
        ffInfo.resize(3);
        for(const auto &current : leptonCurrent) {
#ifdef ENABLE_BSM
            ffInfo[0][current.first] = p_sherpa -> FormFactors(PID::proton(), current.first);
            ffInfo[1][current.first] = p_sherpa -> FormFactors(PID::neutron(), current.first);
            ffInfo[2][current.first] = p_sherpa -> FormFactors(PID::carbon(), current.first);
#else
            // TODO: Define values somewhere
            ffInfo[0][current.first] = SMFormFactor.at({PID::proton(), current.first});
            ffInfo[1][current.first] = SMFormFactor.at({PID::neutron(), current.first});
            ffInfo[2][current.first] = SMFormFactor.at({PID::carbon(), current.first});
#endif
        }
    }
    auto hadronCurrent = m_nuclear -> CalcCurrents(event, ffInfo);

    std::vector<double> amps2(hadronCurrent.size());
    const size_t nlep_spins = leptonCurrent.begin()->second.size();
    const size_t nhad_spins = m_nuclear -> NSpins();
    for(size_t i = 0; i  < nlep_spins; ++i) {
        for(size_t j = 0; j < nhad_spins; ++j) {
            double sign = 1.0;
            std::vector<std::complex<double>> amps(hadronCurrent.size());
            for(size_t mu = 0; mu < 4; ++mu) {
                for(const auto &lcurrent : leptonCurrent) {
                    auto boson = lcurrent.first;
                    for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                        if(hadronCurrent[k].find(boson) != hadronCurrent[k].end()) {
                            amps[k] += sign*lcurrent.second[i][mu]*hadronCurrent[k][boson][j][mu];
                            for(size_t nu = 0; nu < 4; ++nu) {
                                // Calculate W^{\mu\nu}
                                // Calculate L^{h}_{\mu\nu}
                            }
                        }
                    }
                }
                sign = -1.0;
            }
            for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                amps2[k] += std::norm(amps[k]);
            }
        }
    }

    auto lept_in = event.Momentum()[1];
    auto lept_out = event.Momentum().back();
    auto energy_in = lept_in.E();
    auto mass_out = lept_out.M();
    // auto mom_in = lept_in.Vec3().Magnitude();
    auto direction_in = lept_in.Vec3().Unit();
    auto direction_out = lept_out.Vec3().Unit();

    auto hL = FourVector((1/ mass_out) * energy_in * direction_out, (1/ mass_out) * lept_out.Vec3().Magnitude());
    auto hT = FourVector(direction_out.Cross(direction_out.Cross(direction_in)), 0);

    std::map<std::pair<PID, PID>, std::vector<std::array<std::array<std::complex<double>,4>,4>>> hadronTensor;
    std::map<std::pair<PID, PID>, std::array<std::array<std::complex<double>,4>,4>> leptonTensor;
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            for(const auto &lcurrent1 : leptonCurrent) {
                auto boson1 = lcurrent1.first;
                for(const auto &lcurrent2 : leptonCurrent) {
                    auto boson2 = lcurrent2.first;
                    hadronTensor[{boson1, boson2}].resize(hadronCurrent.size());
                    for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                        for(size_t j = 0; j < nhad_spins; ++j) {
                            // Calculate W^{\mu\nu}
                            hadronTensor[{boson1,boson2}][k][mu][nu] += hadronCurrent[k][boson1][j][mu] * std::conj(hadronCurrent[k][boson2][j][nu]);
                        }
                    }
                    for(size_t i = 0; i  < nlep_spins; ++i) {
                        // calculate L^{\mu\nu}
                        leptonTensor[{boson1,boson2}][mu][nu] += (lcurrent1.second[i][mu] * std::conj(lcurrent2.second[i][nu]));
                        //spdlog::info(leptonTensor[{boson1,boson2}][mu][nu]);
                        //spdlog::info(hadronTensor[{boson1,boson2}][0][mu][nu]);
                        //throw;
                    }  
                } 
            }
        }
    }  


    std::array<std::complex<double>, 2> wl{};

    for(const auto &ltensor : leptonTensor) {
        auto bosons = ltensor.first;
        for(size_t k = 0; k < hadronCurrent.size(); ++k) {
            for(size_t mu = 0; mu < 4; ++mu) {
                for(size_t nu = 0; nu < 4; ++nu) {
                    double sign = 1;
                    if ( ((mu == 0) && (nu != 0)) || ((mu != 0) && (nu == 0)) ) {
                        sign = -1;
                    }
                    wl[k] += (leptonTensor[bosons][mu][nu] * sign) * hadronTensor[bosons][k][mu][nu];
                }
            }
        }
    }

    // Contract W^{\mu\nu} with L_{\mu\nu} => Numerator of P_(L,T)
    // getting values for the plot that says theta=16, range of pm .1

    std::vector<double> pL;
    std::vector<double> pT;
 
    double sign = 1;
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            for(size_t alpha = 0; alpha < 4; ++alpha) {
                for(size_t beta = 0; beta < 4; ++beta) {
                    std::complex<double> i = {0,1};
                    pL.resize(hadronCurrent.size());
                    pT.resize(hadronCurrent.size());
                    for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                        if ( ((mu == 0) && (nu != 0)) || ((mu != 0) && (nu == 0)) ) {
                            sign = -1;
                        }
                        if (amps2[1] == 0) {
                            pL[k] = 0;
                            pT[k] = 0;
                        }
                        else {
                            double coupl2 = pow(Constant::ee/(Constant::sw*sqrt(2)), 2);
                            double prefact = coupl2/pow(Constant::MW, 4);
                            pL[k] += prefact * real((mass_out*(hL[mu]*lept_in[nu] + lept_in[mu]*hL[nu] - sign*lept_in*hL 
                            - LeviCivita(mu,nu,alpha,beta)*hL[alpha]*lept_in[beta]*i)*hadronTensor[{-24,-24}][k][mu][nu]))/(real(wl[k]));
                            pT[k] += prefact * real((mass_out*(hT[mu]*lept_in[nu] + lept_in[mu]*hT[nu] - sign*lept_in*hT 
                            + LeviCivita(mu,nu,alpha,beta)*hT[alpha]*lept_in[beta]*i)*hadronTensor[{-24,-24}][k][mu][nu]))/(real(wl[k]));
                        }   
                    }         
                }                   
            }
        }
    }

    event.PolarizationL() = pL;
    event.PolarizationT() = pT; 


    //spdlog::info("pL then pT");
    //spdlog::info(pL[1]);
    //spdlog::info(pT[1]);
    


    double spin_avg = 1;
    if(!ParticleInfo(m_leptonicProcess.m_ids[0]).IsNeutrino()) spin_avg *= 2;
    if(m_nuclear -> NSpins() > 1) spin_avg *= 2;

    // TODO: Correct this flux
    // double flux = 4*sqrt(pow(event.Momentum()[0]*event.Momentum()[1], 2) 
    //                      - event.Momentum()[0].M2()*event.Momentum()[1].M2());
    double mass = ParticleInfo(m_leptonicProcess.m_states.begin()->first[0]).Mass();
    double flux = 2*event.Momentum()[1].E()*2*sqrt(event.Momentum()[0].P2() + mass*mass);
    static constexpr double to_nb = 1e6;
    std::vector<double> xsecs(hadronCurrent.size());
    for(size_t i = 0; i < hadronCurrent.size(); ++i) {
        xsecs[i] = amps2[i]*Constant::HBARC2/spin_avg/flux*to_nb;
        spdlog::debug("Xsec[{}] = {}", i, xsecs[i]);
    }

    return xsecs;
}

bool HardScattering::FillEvent(Event& event, const std::vector<double> &xsecs) const {
    if(!m_nuclear -> FillNucleus(event, xsecs))
        return false;
   
    event.InitializeLeptons(m_leptonicProcess);
    event.InitializeHadrons(m_leptonicProcess);

    return true;
}
