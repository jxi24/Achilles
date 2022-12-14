#include "catch2/catch.hpp" 

#include "Achilles/Constants.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/FourVector.hh"

TEST_CASE("HadronicMapper", "[PhaseSpace]") {
    SECTION("Forward Map") {
        SECTION("Coherent") {
            auto mapper = achilles::CoherentMapper::Construct(0);
            std::vector<achilles::FourVector> mom(1);
            std::vector<double> ran(mapper -> NDims());
            mapper -> GeneratePoint(mom, ran);
            std::vector<double> ran2(mapper -> NDims());
            auto wgt = mapper -> GenerateWeight(mom, ran2);
            CHECK(wgt == 1);
            CHECK(ran == ran2);
        }
        SECTION("QESpectral") {
            auto mapper = achilles::QESpectralMapper::Construct(0);
            std::vector<achilles::FourVector> mom = {{}, {1000, 0, 0, 1000}};
            std::vector<double> ran = {0.5, 0.5, 0.5, 0.5};
            mapper -> SetMasses({0, 0, 0, 0});
            mapper -> GeneratePoint(mom, ran);
            std::vector<double> ran2(4);
            mapper -> GenerateWeight(mom, ran2);
            CHECK(ran == ran2);
            // TODO: Validate wgt
        }
    }

    SECTION("Reverse Map") {
        SECTION("Coherent") {
            auto mapper = achilles::CoherentMapper::Construct(0);
            std::vector<achilles::FourVector> mom = {{achilles::ParticleInfo(achilles::PID::carbon()).Mass(), 0, 0, 0}};
            std::vector<double> ran(mapper -> NDims());
            mapper -> GenerateWeight(mom, ran);
            std::vector<achilles::FourVector> mom2(1);
            mapper -> GeneratePoint(mom2, ran);
            CHECK(mom == mom2);
        }
        SECTION("QESpectral") {
            auto mapper = achilles::QESpectralMapper::Construct(0);
            std::vector<achilles::FourVector> mom = {{achilles::Constant::mN - 20, 400*sin(M_PI/4)*cos(M_PI/6),
                                                      400*sin(M_PI/4)*sin(M_PI/6), 400*cos(M_PI/4)},
                                                     {1000, 0, 0, 1000}};
            std::vector<double> ran(4);
            mapper -> SetMasses({0, 0, 0, 0});
            mapper -> GenerateWeight(mom, ran);
            std::vector<achilles::FourVector> mom2(2);
            mom2[1] = mom[1];
            mapper -> GeneratePoint(mom2, ran);
            CHECK(mom[0].Px() == Approx(mom2[0].Px()));
            CHECK(mom[0].Py() == Approx(mom2[0].Py()));
            CHECK(mom[0].Pz() == Approx(mom2[0].Pz()));
            CHECK(mom[0].E() == Approx(mom2[0].E()));
        }
    }
}
