#include <iostream>
#include <utility>

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

    // basic properties calculation
    auto lept_in = event.Momentum()[1];
    auto lept_out = event.Momentum().back();
    auto energy_in = lept_in.E();
    auto energy_out = lept_out.E();
    auto mom_in = lept_in.Vec3().Magnitude();
    auto mom_out = lept_out.Vec3().Magnitude();
    auto direction_in = lept_in.Vec3().Unit();
    auto direction_out = lept_out.Vec3().Unit();
    auto mass_in = lept_in.M();
    auto mass_out = lept_out.M();

    // print + check basic properties
    /* spdlog::info("{}", lept_in);
    spdlog::info("{}", lept_out);
    spdlog::info("{}", energy_in);
    spdlog::info("{}", energy_out);
    spdlog::info("{}", mom_in);
    spdlog::info("{}", mom_out);
    spdlog::info("{}", direction_in);
    spdlog::info("{}", direction_out);
    spdlog::info("{}", mass_in);
    spdlog::info("{}", mass_out);
    throw; */

    // h_T calculation
    auto ht_vec = direction_out.Cross(direction_out.Cross(direction_in));
    auto ht = FourVector(ht_vec, 0);

    // h_L calculation
    auto hl_vec = energy_in * direction_out;
    auto hl = FourVector(hl_vec, mom_out) / mass_out;

    auto hadronCurrent = m_nuclear -> CalcCurrents(event, ffInfo);

    // print + check h
    /* spdlog::info("{}", hl);
    spdlog::info("{}", ht);
    throw; */

    // amps2[k] calculations
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
                        }
                    }
                }
                sign = -1.0;
            }
            for(size_t k = 0; k < hadronCurrent.size(); ++k)
                amps2[k] += std::norm(amps[k]);
        }
    }

    // calculate spin density matrix: amplitude (depends on k, i, j - nucleons and spins of all particles)

    // W and L calculations
    std::map<std::pair<PID, PID>, std::vector<std::array<std::array<std::complex<double>,4>,4>>> hadronTensor;
    std::map<std::pair<PID, PID>, std::array<std::array<std::complex<double>,4>,4>> leptonTensor;

    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            // loops for gauge bosons
            for(const auto &lcurrent1 : leptonCurrent) {
                auto boson1 = lcurrent1.first;
                for(const auto &lcurrent2 : leptonCurrent) {
                    auto boson2 = lcurrent2.first;
                    hadronTensor[ {boson1, boson2}] .resize(hadronCurrent.size());
                    // loop for nucleon
                    for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                        if( (hadronCurrent[k].find(boson1) != hadronCurrent[k].end()) &&
                            (hadronCurrent[k].find(boson2) != hadronCurrent[k].end()) ) {
                            // loop for spin - hadronic side
                            for(size_t j = 0; j < nhad_spins; ++j) {
                                // Calculate W^{\mu\nu}
                                hadronTensor[ {boson1, boson2} ][k][mu][nu] += hadronCurrent[k][boson1][j][mu] * std::conj(hadronCurrent[k][boson2][j][nu]);
                                // print + check hadronCurrent components
                                /* spdlog::info("{}", hadronCurrent[k][boson1][j][mu]);
                                spdlog::info("{}", std::conj(hadronCurrent[k][boson2][j][nu])); */
                            }
                        }
                    }
                    // loop for spins - leptonic side
                    for(size_t i = 0; i  < nlep_spins; ++i) {
                        // Calculate L^{h}_{\mu\nu}
                        leptonTensor[ {boson1, boson2} ][mu][nu] += lcurrent1.second[i][mu] * std::conj(lcurrent2.second[i][nu]); 
                    }
                }   
            }
        }    
    } 

    // print + check to see if W and amps are 0
    /* for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            for (size_t k = 0; k < hadronCurrent.size(); ++k) {
                if ( (amps2[k] != 0) && (amps2[k] == amps2[k]) ) {
                    spdlog::info("{}", hadronTensor[ {-24,-24} ][k][mu][nu]);
                }
            }
        }
    }

    for (size_t k = 0; k < hadronCurrent.size(); ++k) {
        spdlog::info("{}", amps2[k]);
    }
    throw; */
        

    // compare to eqn 2
    // check to see if you got the right answer: contract hadronic tensor with leptonic tensor
    // another set of for loops over k, mu and nu
    // multiply L_munu by W_kmunu
    // make sure appropriately handle the sign
    // mu=nu=0 +1
    // mu=0 or nu=0 -1
    // every other case is +1
    // L^alphabeta * W^munu * g_alphamu * g_betanu
    // for each k should reproduce amps2[k]

    // amps is a number for protons and a number for neutrons
    // for taus, one is zero

    // 4 combos of spins
    // L_h by summing the 2 +'s with a + and the 2 -'s with a -

    // contraction calculations to check tensor calculations

    std::array<std::complex<double>, 2> contraction;

    for(const auto &ltensor : leptonTensor) {
        auto bosons = ltensor.first;
        for(size_t k = 0; k < hadronCurrent.size(); ++k) {
            for(size_t mu = 0; mu < 4; ++mu) {
                for(size_t nu = 0; nu < 4; ++nu) {
                    double mult = 1;
                    if ( ((mu == 0) && (nu != 0)) || ((mu != 0) && (nu == 0)) ) {
                        mult = -1;
                    }
                    contraction[k] += mult * leptonTensor[bosons][mu][nu] * hadronTensor[bosons][k][mu][nu];
                }
            }
        }
    }

    // print + check if contractions match amps2[k]
    /* if ( (amps2[0] != 0) && (amps2[0] == amps2[0]) ) {
        // spdlog::info("{}", k);
        spdlog::info("{}", amps2[0]);
        spdlog::info("{}", contraction[0]);
    } */

    // print + check amps2[k] if valid
    /* if ( (amps2[0] != 0) && (amps2[0] == amps2[0]) ) {
        spdlog::info("{}", "Denominator:");
        spdlog::info("{}", amps2[0]);
    } */
    // throw;

    // contractions and amps2[k] now agree!

    // Equation 20 calculatios
    // Contract W^{\mu\nu} with L_{\mu\nu} => Numerator of P_(L,T)

    // calculate hk
    std::vector<std::array<std::array<std::complex<double>,4>,4>> hkkh(2);
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            // calculate hk_l + kh_l
            hkkh[0][mu][nu] = hl[mu] * lept_in[nu] + lept_in[mu] * hl[nu];
            // calculate hk_t + kh_t
            hkkh[1][mu][nu] = ht[mu] * lept_in[nu] + lept_in[mu] * ht[nu];
        }
    }
    
    // calculate gkh
    std::vector<std::array<std::array<std::complex<double>,4>,4>> gkh(2);
    for(size_t mu = 0; mu < 4; ++mu) {
        double g = -1;
        if (mu == 0) {
            g = 1;
        }
        // calculate gkh_l
        gkh[0][mu][mu] = g * lept_in * hl;
        // calculate gkh_t
        gkh[1][mu][mu] = g * lept_in * ht;
    }

    // check other values of gkh initiated to zero
    /* for(size_t mu = 0; mu < 4; ++mu) {
        for (size_t nu = 0; nu < 4; ++nu) {
            spdlog::info("{}", gkh[0][mu][nu]);
        }
    }
    throw; */

    // calculate iehk
    std::vector<std::array<std::array<std::complex<double>,4>,4>> iehk(2);
    std::complex<double> i = (0, 1);
    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            for(size_t alpha = 0; alpha < 4; ++alpha) {
                for(size_t beta = 0; beta < 4; ++beta) {
                    // calculate ieh_l * k
                    iehk[0][mu][nu] += i * (std::complex<double>)(LeviCivita(mu, nu, alpha, beta)) * hl[alpha] * lept_in[beta];
                    // print + check iehk components
                    /* if ( (amps2[1] != 0) && (amps2[1] == amps2[1]) ) {
                        spdlog::info("{}", (std::complex<double>)(LeviCivita(mu, nu, alpha, beta)));
                        spdlog::info("{}", hl[alpha]);
                        spdlog::info("{}", ht[alpha]);
                        spdlog::info("{}", lept_in[beta]);
                        spdlog::info("{}", amps2[1]);
                    } */
                    // calculate ieh_t * k
                    iehk[1][mu][nu] += i * (std::complex<double>)(LeviCivita(mu, nu, alpha, beta)) * ht[alpha] * lept_in[beta];
                }
            }
        }
    }
    // throw;

    // calculate numerator
    std::vector<std::array<std::complex<double>, 2>> p_num(2);

    // calculate prefactor and coupling
    double coupl2 = pow(Constant::ee/(Constant::sw*sqrt(2)), 2);
    double prefact = coupl2/pow(Constant::MW, 4);

    for(size_t mu = 0; mu < 4; ++mu) {
        for(size_t nu = 0; nu < 4; ++nu) {
            double mult = 1;
            if ( ((mu == 0) && (nu != 0)) || ((mu != 0) && (nu == 0)) ) {
                mult = -1;
            }
            for(size_t k = 0; k < hadronCurrent.size(); ++k) {
                // newest version - SEGFAULT BEING CAUSED HERE BY ACCESSING HADRON TENSOR
                // print + check values
                /* spdlog::info("{}", mass_out);
                spdlog::info("{}", hkkh[0][mu][nu]);
                spdlog::info("{}", gkh[0][mu][nu]);
                spdlog::info("{}", iehk[0][mu][nu]);
                spdlog::info("{}", hadronTensor[{-24, -24}][k][mu][nu]);
                spdlog::info("{}", mass_out);
                spdlog::info("{}", hkkh[1][mu][nu]);
                spdlog::info("{}", gkh[1][mu][nu]);
                spdlog::info("{}", iehk[1][mu][nu]);
                spdlog::info("{}", hadronTensor[{-24, -24}][k][mu][nu]); */
                if ( (amps2[k] != 0) && (amps2[k] == amps2[k]) ) {
                    p_num[0][k] += mult * prefact * mass_out * (hkkh[0][mu][nu] - gkh[0][mu][nu] + iehk[0][mu][nu]) * hadronTensor[{24, 24}][k][mu][nu];
                    p_num[1][k] += mult * prefact * mass_out * (hkkh[1][mu][nu] - gkh[1][mu][nu] - iehk[1][mu][nu]) * hadronTensor[{24, 24}][k][mu][nu];
                }
                else {
                    p_num[0][k] = 0;
                    p_num[1][k] = 0;
                }
                
            }
            // old version no longer needed because we only need to consider k = 1
            /* for(size_t k = 0; k < hadronCurrent.size(); ++k){
                // calculate numerator for p_l
                p_num[0] += mass_out * (hkkh[0][mu][nu] - gkh[0][mu][nu] + iehk[0][mu][nu]) * hadronTensor[{-24, -24}][k][mu][nu];
                // print + check equation components
                /* spdlog::info("{}", "Numerator Terms:");
                spdlog::info("{}", mass_out);
                spdlog::info("{}", hkkh[0][mu][nu]);
                spdlog::info("{}", gkh[0][mu][nu]);
                spdlog::info("{}", iehk[0][mu][nu]);
                spdlog::info("{}", hadronTensor[{-24, -24}][k][mu][nu]);
                // calculate numerator for p_k
                p_num[1] += mass_out * (hkkh[1][mu][nu] - gkh[1][mu][nu] - iehk[1][mu][nu]) * hadronTensor[{-24, -24}][k][mu][nu];
                // print + check equation components
                spdlog::info("{}", mass_out);
                spdlog::info("{}", hkkh[1][mu][nu]);
                spdlog::info("{}", gkh[1][mu][nu]);
                spdlog::info("{}", iehk[1][mu][nu]);
                spdlog::info("{}", hadronTensor[{-24, -24}][k][mu][nu]);
            }  */
            // old version replaced because need to consider k = 0
            /* p_num[0] += mass_out * (hkkh[0][mu][nu] - gkh[0][mu][nu] + iehk[0][mu][nu]) * hadronTensor[{-24, -24}][1][mu][nu];
            p_num[1] += mass_out * (hkkh[1][mu][nu] - gkh[1][mu][nu] - iehk[1][mu][nu]) * hadronTensor[{-24, -24}][1][mu][nu];
            // print + check equation components
            if ( (amps2[1] != 0) && (amps2[1] == amps2[1]) ) {
                spdlog::info("{}", "Numerator Terms:");
                spdlog::info("{}", mass_out);
                spdlog::info("{}", hkkh[0][mu][nu]);
                spdlog::info("{}", gkh[0][mu][nu]);
                spdlog::info("{}", iehk[0][mu][nu]);
                spdlog::info("{}", hadronTensor[{-24, -24}][1][mu][nu]);
                spdlog::info("{}", hkkh[1][mu][nu]);
                spdlog::info("{}", gkh[1][mu][nu]);
                spdlog::info("{}", iehk[1][mu][nu]); 
            } */
            // print + check W
            /* if ( (amps2[1] != 0) && (amps2[1] == amps2[1]) ) {
                spdlog::info("{}", hadronTensor[{-24, -24}][1][mu][nu]);
            } */
        }
    }

    // print + check num
    /* if ( (amps2[1] != 0) && (amps2[1] == amps2[1]) ) {
        spdlog::info("{}", "Numerators:");
        spdlog::info("{}", p_num[0]);
        spdlog::info("{}", p_num[1]);
    } */

    // calculate, print + check final polarization vector
    std::vector<std::array<std::complex<double>,2>> p(hadronCurrent.size());
    // spdlog::info("{}", "Results:");
    // loop over polarization
    for(size_t i = 0; i < 2; ++i) {
        for(size_t k = 0; k < hadronCurrent.size(); ++k) {
            if ( (amps2[k] != 0) && (amps2[k] == amps2[k]) ) {
                p[i][k] = p_num[i][k] / amps2[k];
                // spdlog::info("{}", p[i][k]);
                // spdlog::info("{}", "amps2[k] is valid");
            }
            else if (amps2[k] != amps2[k]) {
                p[i][k] = 0;
                // spdlog::info("{}", "amps2[k] is nan");
                // throw;
            }
            else {
                p[i][k] = 0;
                // spdlog::info("{}", p[i][k]);
                // spdlog::info("{}", "amps2[k] is nan");
                // throw;
            }
        }  
    } 

    // print + check final polarization vector
    // p_l
    /* spdlog::info("{}", "p[0][0]=");
    spdlog::info("{}", p[0][0]);
    spdlog::info("{}", "p[0][1]=");
    spdlog::info("{}", p[0][1]);
    // p_t
    spdlog::info("{}", "p[1][0]=");
    spdlog::info("{}", p[1][0]);
    spdlog::info("{}", "p[1][1]=");
    spdlog::info("{}", p[1][1]);
    throw; */

    // connect final polarization vector to event code
    event.set_polarization_l({p[0][0].real(), p[0][1].real()});
    event.set_polarization_t({p[1][0].real(), p[1][1].real()});

    // send amps2[k] info to event gen
    event.set_amps2( {amps2[0], amps2[1]} );
    // event gen only exports value if amps2[k] is not nan

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
