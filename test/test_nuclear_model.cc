#include "catch2/catch.hpp"

#include "Achilles/NuclearModel.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"

#include "yaml-cpp/yaml.h"
#include "mock_classes.hh"
#include "Approx.hh"

using achilles::operator""_GeV;

TEST_CASE("CoherentModel", "[NuclearModel]") {
    YAML::Node config = YAML::Load(R"config(
Nucleus:
  Name: 12C
  Binding: 0
  Fermi Momentum: 0
  FermiGas: Local
  Density:
    File: data/c12.prova.txt
)config");
    YAML::Node ff = YAML::Load("vector: dummy\naxial: dummy\ncoherent: dummy\ndummy: dummy2");

    auto form_factor = std::make_unique<MockFormFactor>();

    MockFormFactorBuilder builder;
    trompeloeil::sequence seq;
    REQUIRE_CALL(builder, Vector(ff["vector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, AxialVector(ff["axial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, Coherent(ff["coherent"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));

    SECTION("Allowed States are valid") {
        // Require to build here or else form_factor is moved before expectations are set in next test
        REQUIRE_CALL(builder, build())
            .TIMES(1)
            .IN_SEQUENCE(seq)
            .LR_RETURN(std::move(form_factor));
        achilles::Coherent model(config, ff, builder);

        achilles::Process_Info info;
        info.m_ids = {achilles::PID::electron(), achilles::PID::electron()};
        model.AllowedStates(info);
        CHECK(info.m_states[{achilles::PID::carbon()}] == std::vector<achilles::PID>{achilles::PID::carbon()});

        achilles::Process_Info invalid;
        invalid.m_ids = {achilles::PID::electron(), achilles::PID::nu_electron()};
        CHECK_THROWS_WITH(model.AllowedStates(invalid), "Achilles::NuclearError:Coherent: Requires charge 0, but found charge 1");
    }

    SECTION("CalcCurrents") {
        std::vector<achilles::FourVector> momentum = {{11.178_GeV, 0, 0, 0},
                                                      {4.159051495317648_GeV, 0, 0, 4.159051495317648_GeV},
                                                      {11.179382065495107_GeV, 0.11140855999009017_GeV,
                                                       0.1334001474175966_GeV, 0.0263039872165316_GeV},
                                                      {3.195267141839163_GeV, -0.19106532624117684_GeV,
                                                       0.036750024067231254_GeV, 3.189337797185557_GeV},
                                                      {0.8015109318750572_GeV, 0.07305947032271103_GeV,
                                                       -0.13285393143782712_GeV, 0.7870399739612345_GeV},
                                                      {0.16089135610832217_GeV, 0.006597295928375643_GeV,
                                                       -0.03729624004700072_GeV, 0.15636973695432532_GeV}};
        MockEvent e;
        const MockEvent& event = e;
        REQUIRE_CALL(event, Momentum())
            .TIMES(AT_LEAST(5))
            .LR_RETURN((momentum));

        double Q2 = -(momentum[1] - momentum[3] - momentum[4] - momentum[5]).M2();
        achilles::FormFactor::Values value;
        value.Fcoh = 1;
        REQUIRE_CALL(*form_factor, call_op(Q2))
            .TIMES(1)
            .LR_RETURN((value));

        // Require to build here or else form_factor is moved before expectations are set
        REQUIRE_CALL(builder, build())
            .TIMES(1)
            .IN_SEQUENCE(seq)
            .LR_RETURN(std::move(form_factor));
        achilles::Coherent model(config, ff, builder);

        std::vector<achilles::NuclearModel::FFInfoMap> info_map(3);
        info_map[2][achilles::PID::carbon()] = {achilles::FormFactorInfo{achilles::FormFactorInfo::Type::FCoh, 1}};
        auto results = model.CalcCurrents(event, info_map);
        std::vector<std::complex<double>> expected = {momentum[0][0]+momentum[2][0],
                                                      momentum[0][1]+momentum[2][1],
                                                      momentum[0][2]+momentum[2][2],
                                                      momentum[0][3]+momentum[2][3]};
        CHECK(results[0][achilles::PID::carbon()][0] == expected);
    }

    SECTION("Properly fill event") {
        // Require to build here or else form_factor is moved before expectations are set
        REQUIRE_CALL(builder, build())
            .TIMES(1)
            .IN_SEQUENCE(seq)
            .LR_RETURN(std::move(form_factor));
        achilles::Coherent model(config, ff, builder);

        std::vector<achilles::Particle> nucleons{};
        auto nucleus = std::make_shared<MockNucleus>();
        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(3)
            .LR_RETURN((nucleons));

        std::vector<achilles::FourVector> momentum = {{11.178_GeV, 0, 0, 0},
                                                      {4.159051495317648_GeV, 0, 0, 4.159051495317648_GeV},
                                                      {11.179382065495107_GeV, 0.11140855999009017_GeV,
                                                       0.1334001474175966_GeV, 0.0263039872165316_GeV},
                                                      {3.195267141839163_GeV, -0.19106532624117684_GeV,
                                                       0.036750024067231254_GeV, 3.189337797185557_GeV},
                                                      {0.8015109318750572_GeV, 0.07305947032271103_GeV,
                                                       -0.13285393143782712_GeV, 0.7870399739612345_GeV},
                                                      {0.16089135610832217_GeV, 0.006597295928375643_GeV,
                                                       -0.03729624004700072_GeV, 0.15636973695432532_GeV}};

        achilles::Event event;
        event.CurrentNucleus() = nucleus;
        event.Momentum() = momentum;

        std::vector<double> xsecs = {10};
        auto pass = model.FillNucleus(event, xsecs);
        CHECK(pass == true);
    }
}

TEST_CASE("QESpectralModel", "[NuclearModel]") {
    YAML::Node config = YAML::Load(R"config(
NuclearModel:
  SpectralP: data/pke12_tot.data
  SpectralN: data/pke12_tot.data
  Ward: false
)config");
    YAML::Node ff = YAML::Load("vector: dummy\naxial: dummy\ncoherent: dummy\ndummy: dummy2");

    auto form_factor = std::make_unique<MockFormFactor>();

    MockFormFactorBuilder builder;
    trompeloeil::sequence seq;
    REQUIRE_CALL(builder, Vector(ff["vector"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, AxialVector(ff["axial"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));
    REQUIRE_CALL(builder, Coherent(ff["coherent"].as<std::string>(), ff["dummy"]))
        .TIMES(1)
        .IN_SEQUENCE(seq)
        .LR_RETURN(std::ref(builder));

    SECTION("Allowed States are valid") {
        // Require to build here or else form_factor is moved before expectations are set in next test
        REQUIRE_CALL(builder, build())
            .TIMES(1)
            .IN_SEQUENCE(seq)
            .LR_RETURN(std::move(form_factor));
        achilles::QESpectral model(config, ff, builder);

        SECTION("Neutral Current") {
            achilles::Process_Info info;
            info.m_ids = {achilles::PID::electron(), achilles::PID::electron()};
            model.AllowedStates(info);
            CHECK(info.m_states[{achilles::PID::proton()}] == std::vector<achilles::PID>{achilles::PID::proton()});
            CHECK(info.m_states[{achilles::PID::neutron()}] == std::vector<achilles::PID>{achilles::PID::neutron()});
        }

        SECTION("-1 Nuclear Charge") {
            achilles::Process_Info info;
            info.m_ids = {-achilles::PID::nu_electron(), -achilles::PID::electron()};
            model.AllowedStates(info);
            CHECK(info.m_states[{achilles::PID::proton()}] == std::vector<achilles::PID>{achilles::PID::neutron()});
        }

        SECTION("+1 Nuclear Charge") {
            achilles::Process_Info info;
            info.m_ids = {achilles::PID::nu_electron(), achilles::PID::electron()};
            model.AllowedStates(info);
            CHECK(info.m_states[{achilles::PID::neutron()}] == std::vector<achilles::PID>{achilles::PID::proton()});
        }

        SECTION("-2 Nuclear Charge") {
            achilles::Process_Info invalid;
            invalid.m_ids = {achilles::PID::electron(), -achilles::PID::electron()};
            CHECK_THROWS_WITH(model.AllowedStates(invalid), "Achilles::NuclearError:Quasielastic: Requires |charge| < 2, but found |charge| 2");
        }

        SECTION("+2 Nuclear Charge") {
            achilles::Process_Info invalid;
            invalid.m_ids = {-achilles::PID::electron(), achilles::PID::electron()};
            CHECK_THROWS_WITH(model.AllowedStates(invalid), "Achilles::NuclearError:Quasielastic: Requires |charge| < 2, but found |charge| 2");
        }
    }

    SECTION("CalcCurrents") {
        std::vector<achilles::FourVector> momentum = {{0.8334268643628409_GeV, 0.0841386014098756_GeV,
                                                       0.35434104526508325_GeV, -0.25207280069997196_GeV},
                                                      {4_GeV, 0, 0, 4_GeV},
                                                      {3.9936896971024103_GeV, 0.8919304531816128_GeV,
                                                       0.5832575462220676_GeV, 3.732825216669567_GeV},
                                                      {0.8397371672604309_GeV, -0.8077918517717372_GeV,
                                                       -0.22891650095698435_GeV, 0.01510198263046103_GeV}};
        MockEvent e;
        const MockEvent& event = e;
        REQUIRE_CALL(event, Momentum())
            .TIMES(AT_LEAST(5))
            .LR_RETURN((momentum));

        double Q2 = -(momentum[1] - momentum[3]).M2()/1.0_GeV/1.0_GeV;
        achilles::FormFactor::Values value;
        value.F1p = 1;
        value.F1n = 1;
        value.F2p = 1;
        value.F2n = 1;
        value.FA = 1;
        REQUIRE_CALL(*form_factor, call_op(Q2))
            .TIMES(1)
            .LR_RETURN((value));

        // Require to build here or else form_factor is moved before expectations are set in next test
        REQUIRE_CALL(builder, build())
            .TIMES(1)
            .IN_SEQUENCE(seq)
            .LR_RETURN(std::move(form_factor));
        achilles::QESpectral model(config, ff, builder);

        std::vector<achilles::NuclearModel::FFInfoMap> info_map(3);
        info_map[0][achilles::PID::photon()] = {achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F1p, 1},
                                                achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F2p, 1},
                                                achilles::FormFactorInfo{achilles::FormFactorInfo::Type::FA, 1}};
        info_map[1][achilles::PID::photon()] = {achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F1n, 1},
                                                achilles::FormFactorInfo{achilles::FormFactorInfo::Type::F2n, 1},
                                                achilles::FormFactorInfo{achilles::FormFactorInfo::Type::FA, 1}};
        auto results = model.CalcCurrents(event, info_map);
        std::vector<std::vector<std::complex<double>>> expected = {{{0.00018993,0.0012343},
                                                                    {-0.00170272,4.6931e-05},
                                                                    {0.000416924,-0.00136515},
                                                                    {0.000520202,0.00136208}},
                                                                   {{-0.0043949,0.000799289},
                                                                    {-0.00120857,0.000629831},
                                                                    {-0.0015762,0.000136311},
                                                                    {-0.00393553,0.000641262}},
                                                                   {{0.00200181,0.000364063},
                                                                    {0.000714529,0.0027743},
                                                                    {0.00336446,-0.000357726},
                                                                    {0.000134006,-0.000455973}},
                                                                   {{-0.000253546,0.00164772},
                                                                    {-0.00670183,0.00152014},
                                                                    {0.0013261,0.00679448},
                                                                    {0.00103612,0.000901872}},
                                                                   {{0.00018993,0.0012343},
                                                                    {-0.00170272,4.6931e-05},
                                                                    {0.000416924,-0.00136515},
                                                                    {0.000520202,0.00136208}},
                                                                   {{-0.0043949,0.000799289},
                                                                    {-0.00120857,0.000629831},
                                                                    {-0.0015762,0.000136311},
                                                                    {-0.00393553,0.000641262}},
                                                                   {{0.00200181,0.000364063},
                                                                    {0.000714529,0.0027743},
                                                                    {0.00336446,-0.000357726},
                                                                    {0.000134006,-0.000455973}},
                                                                   {{-0.000253546,0.00164772},
                                                                    {-0.00670183,0.00152014},
                                                                    {0.0013261,0.00679448},
                                                                    {0.00103612,0.000901872}}};
        for(size_t i = 0; i < 4; ++i) {
            REQUIRE_THAT(results[0][achilles::PID::photon()][i], VectorComplexApprox(expected[i]).margin(1e-5));
        }
    }

    SECTION("Properly fill the event") {
        // Require to build here or else form_factor is moved before expectations are set in next test
        REQUIRE_CALL(builder, build())
            .TIMES(1)
            .IN_SEQUENCE(seq)
            .LR_RETURN(std::move(form_factor));
        achilles::QESpectral model(config, ff, builder);

        std::vector<achilles::Particle> nucleons{achilles::PID::proton(), achilles::PID::neutron()};
        auto nucleus = std::make_shared<MockNucleus>();
        REQUIRE_CALL(*nucleus, Nucleons())
            .TIMES(5) // Should be 2*nucleons.size()+1
            .LR_RETURN((nucleons));
        achilles::Event event(1);
        event.MatrixElementWgts().resize(nucleons.size());
        event.CurrentNucleus() = nucleus;

        std::vector<double> xsecs = {10, 20};
        auto pass = model.FillNucleus(event, xsecs);
        CHECK(pass == true);
        event.CalcWeight();
        CHECK(event.Weight() == xsecs[0] + xsecs[1]);
    }
}
