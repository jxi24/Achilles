#ifndef EVENT_HH
#define EVENT_HH

#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "Achilles/Achilles.hh"
#include "Achilles/HardScatteringEnum.hh"
#include "Achilles/NuclearRemnant.hh"
#include "Achilles/ProcessInfo.hh"

namespace achilles {

class PID;
class FourVector;
class Particle;
class Nucleus;
class Beam;
class NuclearModel;

struct InitialState {
    std::shared_ptr<Beam> beam;
    std::shared_ptr<Nucleus> nucleus;
};

using vParticles = std::vector<Particle>;
using vMomentum = std::vector<FourVector>;

class Event {
    public:
        Event(double vWgt = 0) : m_vWgt{vWgt} {}
        Event(std::shared_ptr<Nucleus>, 
              std::vector<FourVector>, double);
        MOCK ~Event() = default;

        void SetHardScatteringType(HardScatteringType type) { m_type = type; }
        MOCK void InitializeLeptons(const Process_Info&);
        MOCK void InitializeHadrons(const Process_Info&);
        void Finalize();

        MOCK const NuclearRemnant &Remnant() const { return m_remnant; }

        MOCK const vMomentum &Momentum() const { return m_mom; }
        MOCK vMomentum &Momentum() { return m_mom; }

        const double &MatrixElementWgt(size_t i) const { return m_me[i]; }
        double &MatrixElementWgt(size_t i) { return m_me[i]; }

        const std::vector<double> &MatrixElementWgts() const { return m_me; }
        std::vector<double> &MatrixElementWgts() { return m_me; }

        bool TotalCrossSection();
        size_t SelectNucleon() const;

        const std::shared_ptr<Nucleus>& CurrentNucleus() const { return m_nuc; }
        MOCK std::shared_ptr<Nucleus>& CurrentNucleus() { return m_nuc; }

        MOCK vParticles Particles() const;
        const vParticles& Hadrons() const;
        MOCK vParticles& Hadrons();
        const vParticles& Leptons() const { return m_leptons; }
        vParticles& Leptons() { return m_leptons; }
        MOCK double Weight() const;
        void SetMEWeight(double wgt) { m_meWgt = wgt; }
        void Rotate(const std::array<double,9>&);

        bool operator==(const Event &other) const {
            return m_type == other.m_type && m_nuc == other.m_nuc
                && m_remnant == other.m_remnant && m_mom == other.m_mom
                && m_me == other.m_me && m_vWgt == other.m_meWgt
                && m_leptons == other.m_leptons;
        }

        // accessor, returns an array of polarization_l for k = 0 and k = 1
        std::array<double, 2> get_polarization_l();
        // accessor, returns an array of polarization_t for k = 0 and k = 1
        std::array<double, 2> get_polarization_t();
        // mutator, modifies an array of polarization_l for k = 0 and k = 1
        void set_polarization_l(std::array<double, 2> new_polarization);
        // mutator, modifies an array of polarization_l for k = 0 and k = 1
        void set_polarization_t(std::array<double, 2> new_polarization);
        // accessor, returns amps2[k]
        std::array<std::complex<double>, 2> get_amps2();
        // mutator, defines amps2[k]
        void set_amps2(std::array<std::complex<double>, 2> new_amps2);

    private:
        std::vector<double> EventProbs() const;

        // bool ValidateEvent(size_t) const;

        HardScatteringType m_type{HardScatteringType::None};
        std::shared_ptr<Nucleus> m_nuc;
        NuclearRemnant m_remnant{};
        vMomentum m_mom{};
        std::vector<double> m_me;
        double m_vWgt{}, m_meWgt{};
        vParticles m_leptons{};
        vParticles m_history{};
        // arrays keeping track of p_l and p_t for k = 0 and k = 1
        std::array<double, 2> m_polarization_l{};
        std::array<double, 2> m_polarization_t{};
        // array keeping track of amps2[k]
        std::array<std::complex<double>, 2> m_amps2{};
};

}

#endif
