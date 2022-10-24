#include "catch2/catch.hpp"

#include "Achilles/CombinedCuts.hh"
#include "Achilles/Cuts.hh"
#include "Achilles/OneParticleCuts.hh"
#include "Achilles/TwoParticleCuts.hh"
#include "Achilles/FourVector.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"

#include "catch_utils.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

#include <iostream>

template<class T>
struct TestCutBase : public achilles::CutBase<T> {
    TestCutBase(YAML::Node node) : achilles::CutBase<T>(node) {}
    bool MakeCut(const T& val) const { return this->CheckCut(val); }
};

TEMPLATE_TEST_CASE("Base Cut", "[Cuts]", double, float, size_t, int) {
    SECTION("Check for valid syntax") {
        YAML::Node node;
        REQUIRE_THROWS_AS(TestCutBase<TestType>(node), std::logic_error);
        node["min"] = 30;
        REQUIRE_NOTHROW(TestCutBase<TestType>(node));
        node["max"] = 100;
        REQUIRE_NOTHROW(TestCutBase<TestType>(node));
        node["range"] = std::vector<TestType>{30, 50};
        REQUIRE_THROWS_AS(TestCutBase<TestType>(node), std::logic_error);
        YAML::Node node2;
        node2["range"] = std::vector<TestType>{30, 50};
        REQUIRE_NOTHROW(TestCutBase<TestType>(node2));
    }

    SECTION("Cuts work as expected") {
        YAML::Node node;
        node["min"] = 10;
        auto cut = TestCutBase<TestType>(node);
        auto val = static_cast<TestType>(GENERATE(take(30, random(0., 100.))));
        CHECK((val > 10) == cut.MakeCut(val)); 
    }
}

TEST_CASE("Single particle cuts", "[Cuts]") {
    auto mom = GENERATE(take(30, randomMomentum(1000)));

    SECTION("Energy Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::EnergyCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.E() > 10));
    }

    SECTION("Momentum Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::MomentumCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.P() > 10));
    }

    SECTION("AngleTheta Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::AngleThetaCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.Theta()*180/M_PI > 10));
    }

    SECTION("Transverse Momentum Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::TransverseMomentumCut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.Pt() > 10));
    }

    SECTION("ETheta2 Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::ETheta2Cut cut(node);
        CHECK(cut.MakeCut(mom) == (mom.E()*pow(mom.Theta(), 2) > 10));
    }
}

TEST_CASE("Two particle cuts", "[Cuts]") {
    auto mom1 = GENERATE(take(30, randomMomentum(1000)));
    auto mom2 = GENERATE(take(30, randomMomentum(1000)));

    SECTION("DeltaTheta Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::DeltaThetaCut cut(node);
        CHECK(cut.MakeCut(mom1, mom2) == (std::acos(mom1.CosAngle(mom2))*180/M_PI > 10));
    }

    SECTION("InvariantMass Cut") {
        YAML::Node node;
        node["min"] = 10;
        achilles::InvariantMassCut cut(node);
        CHECK(cut.MakeCut(mom1, mom2) == ((mom1+mom2).M() > 10));
    }
}

TEST_CASE("CutCollection", "[Cuts]") {
    achilles::CutCollection cuts;
    auto mom1 = GENERATE(take(30, randomMomentum(1000)));
    auto mom2 = GENERATE(take(30, randomMomentum(1000)));
    auto part1 = achilles::Particle(achilles::PID::electron(), mom1);
    auto part2 = achilles::Particle(achilles::PID::muon(), mom2);
    auto part3 = achilles::Particle(achilles::PID::electron(), mom2);
    auto part4 = achilles::Particle(-achilles::PID::electron(), mom2);
    part1.Status() = achilles::ParticleStatus::final_state;
    part2.Status() = achilles::ParticleStatus::final_state;
    part3.Status() = achilles::ParticleStatus::final_state;
    part4.Status() = achilles::ParticleStatus::final_state;

    SECTION("Single PIDs work") {
        YAML::Node node;
        node["min"] = 10;
        auto cut = achilles::CutFactory<achilles::OneParticleCut>::InitializeCut("Energy", node);
        cuts.AddCut({achilles::PID::electron()}, std::move(cut)); 
        CHECK(cuts.EvaluateCuts({part1}) == (mom1.E() > 10));
        CHECK(cuts.EvaluateCuts({part2}) == true);
        CHECK(cuts.EvaluateCuts({part1, part2}) == (mom1.E() > 10));

        auto cut2 = achilles::CutFactory<achilles::TwoParticleCut>::InitializeCut("DeltaTheta", node);
        cuts.AddCut({achilles::PID::electron()}, std::move(cut2));
        bool pass_cuts = (mom1.E() > 10 && mom2.E() > 10 && std::acos(mom1.CosAngle(mom2))*180/M_PI > 10);
        CHECK(cuts.EvaluateCuts({part1, part3}) == pass_cuts);
    }

    SECTION("Sets of PIDs work") {
        YAML::Node node;
        node["min"] = 10;
        auto cut = achilles::CutFactory<achilles::OneParticleCut>::InitializeCut("Energy", node);
        cuts.AddCut({achilles::PID::electron(), achilles::PID::muon()}, std::move(cut)); 
        CHECK(cuts.EvaluateCuts({part1, part2}) == (mom1.E() > 10 && mom2.E() > 10));

        auto cut2 = achilles::CutFactory<achilles::TwoParticleCut>::InitializeCut("DeltaTheta", node);
        cuts.AddCut({achilles::PID::electron(), achilles::PID::muon()}, std::move(cut2));
        bool pass_cuts = (mom1.E() > 10 && mom2.E() > 10 && std::acos(mom1.CosAngle(mom2))*180/M_PI > 10);
        CHECK(cuts.EvaluateCuts({part1, part2}) == pass_cuts);
    }

    SECTION("YAML correctly builds CutCollection") {
        YAML::Node node = YAML::Load(R"node(
        cuts:
            - Type: Energy
              PIDs: [11, -11]
              min: 10
            - Type: InvariantMass
              PIDs: [11, -11]
              range: [80, 100]
            - Type: AngleTheta
              PIDs: 13
              max: 30
            - Type: DeltaTheta
              PIDs: [11, -11]
              range: [[10, 30], [60, 80]]
        )node");

        cuts = node["cuts"].as<achilles::CutCollection>();
        bool pass_cuts = (mom1.E() > 10) && (mom2.E() > 10);
        pass_cuts &= (mom1+mom2).M() > 80 && (mom1+mom2).M() < 100;
        double dtheta = std::acos(mom1.CosAngle(mom2))*180/M_PI;
        pass_cuts &= (dtheta > 10 && dtheta < 30) || (dtheta > 60 && dtheta < 80);
        CHECK(cuts.EvaluateCuts({part1, part4}) == pass_cuts);
    }
}
