#include "Achilles/ProcessInfo.hh"

size_t achilles::Process_Info::Multiplicity() const {
    return m_ids.size() + m_states.begin() -> first.size() + m_states.begin() -> second.size();
}

std::vector<double> achilles::Process_Info::Masses() const {
    std::vector<double> masses;

    // Get Hadronic current masses
    for(size_t i = 0; i < m_states.begin() -> second.size(); ++i) {
        spdlog::trace("Mass[{}] = {}", m_states.begin() -> second[i], ParticleInfo(m_states.begin() -> second[i]).Mass());
        masses.push_back(pow(ParticleInfo(m_states.begin() -> second[i]).Mass(), 2));
    }

    // Get Leptonic current masses
    for(size_t i = 1; i < m_ids.size(); ++i) {
        spdlog::trace("Mass[{}, {}] = {}", m_ids[i], ParticleInfo(m_ids[i]).Name(), ParticleInfo(m_ids[i]).Mass());
        masses.push_back(pow(ParticleInfo(m_ids[i]).Mass(), 2));
    }

    return masses;
}

std::vector<long> achilles::Process_Info::Ids() const {
    std::vector<long> ids;

    // Get initial hadronic id
    ids.push_back(m_states.begin() -> first[0].AsInt());

    // Get inital lepton id
    ids.push_back(m_ids[0].AsInt());

    // Get remaining hadronic ids
    for(size_t i = 0; i < m_states.begin() -> second.size(); ++i)
        ids.push_back(m_states.begin() -> second[i].AsInt());

    // Get remaining leptonic ids
    for(size_t i = 1; i < m_ids.size(); ++i)
        ids.push_back(m_ids[i].AsInt());

    return ids;
}
