#include "Achilles/Event.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Beams.hh"
#include "Achilles/NuclearModel.hh"

using achilles::Event;

Event::Event(std::shared_ptr<Nucleus> nuc,
             std::vector<FourVector> mom, double vwgt)
        : m_nuc{std::move(nuc)}, m_mom{std::move(mom)}, m_vWgt{std::move(vwgt)} {
    m_nuc -> GenerateConfig();
    m_me.resize(m_nuc -> NNucleons());
}

// bool Event::ValidateEvent(size_t imatrix) const {
//     spdlog::trace("Number of momentums = {}", m_mom.size());
//     spdlog::trace("Number of states = {}", m_me[imatrix].inital_state.size() + m_me[imatrix].final_state.size());
//     return m_mom.size() == m_me[imatrix].inital_state.size() + m_me[imatrix].final_state.size();
// }

// TODO: Make this cleaner and handle multiple nucleons in initial state
void Event::InitializeLeptons(const Process_Info &process) {
    // Setup leptons
    bool initial_state = true;
    for(const auto &elm : process.m_mom_map) {
        if(ParticleInfo(elm.second).IsLepton() || ParticleInfo(elm.second).IsVector()) {
            m_leptons.emplace_back(ParticleInfo(elm.second), m_mom[elm.first]);
            if(initial_state) {
                m_leptons.back().Status() = ParticleStatus::initial_state;
                initial_state = false;
            } else {
                m_leptons.back().Status() = ParticleStatus::final_state;
            }
        }
    }
}

void Event::InitializeHadrons(const Process_Info &process) {
    // Coherent scattering is handled inside the coherent class
    if(ParticleInfo(process.m_states.begin()->first[0]).IsNucleus())
        return;

    // Get all hadronic momenta
    std::vector<FourVector> mom;
    for(const auto &elm : process.m_mom_map) {
        if(ParticleInfo(elm.second).IsHadron()) {
            mom.push_back(m_mom[elm.first]);
        }
    }

    // TODO: Update to handle multiple initial and final state particles
    // Initial state setup
    size_t idx = SelectNucleon();
    Particle &initial = m_nuc -> Nucleons()[idx];
    initial.Momentum() = mom.front();
    initial.Status() = ParticleStatus::initial_state;

    // Final state setup
    Particle final(process.m_states.at({initial.ID()})[0], mom.back(),
                   initial.Position(), ParticleStatus::propagating);
    m_nuc -> Nucleons().push_back(final);
}

void Event::Finalize() {
    size_t nA = 0, nZ = 0;
    for(auto it = m_nuc -> Nucleons().begin(); it != m_nuc -> Nucleons().end(); ) {
        if(it -> Status() == ParticleStatus::background) {
            if(it -> ID() == PID::proton()) nZ++;
            nA++;
            it = m_nuc -> Nucleons().erase(it);
        } else {
            ++it;
        }
    }

    m_remnant = NuclearRemnant(nA, nZ);
}

achilles::vParticles Event::Particles() const {
    vParticles result;
    result.insert(result.end(), Hadrons().begin(), Hadrons().end());
    result.insert(result.end(), m_leptons.begin(), m_leptons.end());
    return result;
}

const achilles::vParticles& Event::Hadrons() const {
    return m_nuc -> Nucleons();
}

achilles::vParticles& Event::Hadrons() {
    return m_nuc -> Nucleons();
}

double Event::Weight() const {
    return m_vWgt*m_meWgt;
}

bool Event::TotalCrossSection() {
    m_meWgt = std::accumulate(m_me.begin(), m_me.end(), 0.0, std::plus<>());
    spdlog::debug("Total xsec = {}", m_meWgt);
    return m_meWgt > 0 ? true : false;
}

size_t Event::SelectNucleon() const {
    std::vector<double> probs = EventProbs();
    double rand = Random::Instance().Uniform(0.0, 1.0);
    return static_cast<size_t>(std::distance(probs.begin(),
                               std::lower_bound(probs.begin(), probs.end(), rand)))-1;
}

std::vector<double> Event::EventProbs() const {
    std::vector<double> probs;
    probs.push_back(0.0);
    double cumulative = 0.0;
    for(const auto & m : m_me) {
        cumulative += m;
        probs.emplace_back(cumulative / m_meWgt);
    }
    return probs;
}

void Event::Rotate(const std::array<double,9>& rot_mat) {
    for (auto& particle: m_nuc -> Nucleons()){ particle.Rotate(rot_mat); }
    for (auto& particle: m_leptons){ particle.Rotate(rot_mat); }
}

std::array<double, 2> Event::get_polarization_l() {
    return m_polarization_l;
}

std::array<double, 2> Event::get_polarization_t() {
    return m_polarization_t;
}

void Event::set_polarization_l(std::array<double, 2> new_polarization) {
    m_polarization_l = new_polarization;
}

void Event::set_polarization_t(std::array<double, 2> new_polarization) {
    m_polarization_t = new_polarization;
}

std::array<std::complex<double>, 2> Event::get_amps2() {
    return m_amps2;
}

void Event::set_amps2(std::array<std::complex<double>, 2> new_amps2) {
    m_amps2 = new_amps2;
}


