#ifndef EVENTGEN_HH
#define EVENTGEN_HH

#include "Achilles/CombinedCuts.hh"
#include "Achilles/Histogram.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/QuasielasticTestMapper.hh"
#include "Achilles/Vegas.hh"
#include "Achilles/MultiChannel.hh"

#include <memory>
#include <vector>

namespace YAML {

class Node;

}

namespace achilles {

// Forward declare types
class Event;
class Beam;
class Nucleus;
class Cascade;
class HardScattering;
class EventWriter;

class SherpaMEs;

class EventGen {
    public:
        EventGen(const std::string&, std::vector<std::string>);
        void Initialize();
        // helper function to calculate new error
        double PropogateError(const double p_1, const double p_2, const double error_1, const double error_2) const;
        void GenerateEvents();

    private:
        bool runCascade{false}, outputEvents{false}, doHardCuts{false}, doEventCuts{false};
        bool doRotate{false};
        double GenerateEvent(const std::vector<FourVector>&, const double&);
        bool MakeCuts(Event&);
        // bool MakeEventCuts(Event&);
        void Rotate(Event&);

        std::shared_ptr<Beam> beam;
        std::shared_ptr<Nucleus> nucleus;
        std::shared_ptr<Cascade> cascade;
        std::shared_ptr<HardScattering> scattering;
        CutCollection hard_cuts{};
        // CutCollection event_cuts{};
        MultiChannel integrator;
        Integrand<FourVector> integrand;
        YAML::Node config;

        std::shared_ptr<EventWriter> writer;

        // array keeping track of p_L StatsData objects for k = 0 and k = 1
        std::array<StatsData, 2> Polarization_l;
        // array keeping track of p_T StatsData objects for k = 0 and k = 1
        std::array<StatsData, 2> Polarization_t;
        // array keeping track of amps2[k]
        std::array<double, 2> Amps2;
        // double keeping track of q0
        double Q0;
        // double keeping track of theta
        double Theta;
};

}

#endif
