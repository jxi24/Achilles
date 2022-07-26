#include "Achilles/EventGen.hh"
#include "Achilles/Event.hh"
#include "Achilles/EventWriter.hh"
#include "Achilles/HardScatteringFactory.hh"
#include "Achilles/HardScattering.hh"
#include "Achilles/Logging.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Beams.hh"
#include "Achilles/Cascade.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Units.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/NuclearModel.hh"
#include "Achilles/ComplexFmt.hh"
#include "Achilles/Units.hh"

// TODO: Turn this into a factory to reduce the number of includes
#include "Achilles/PhaseSpaceBuilder.hh"
#include "Achilles/BeamMapper.hh"
#include "Achilles/HadronicMapper.hh"
#include "Achilles/FinalStateMapper.hh"
#include "Achilles/PhaseSpaceMapper.hh"
#include "Achilles/QuasielasticTestMapper.hh"

#include <iostream>
#include <string>
#include <fstream>

#ifdef ENABLE_BSM
#include "plugins/Sherpa/Channels1.hh"
#include "plugins/Sherpa/Channels3.hh"
#include "plugins/Sherpa/SherpaMEs.hh"
#endif

#ifdef ENABLE_HEPMC3
#include "plugins/HepMC3/HepMC3EventWriter.hh"
#endif

#include "yaml-cpp/yaml.h"

achilles::Channel<achilles::FourVector> BuildChannelTest(const YAML::Node &node, std::shared_ptr<achilles::Beam> beam) {
    achilles::Channel<achilles::FourVector> channel;
    channel.mapping = std::make_unique<achilles::QuasielasticTestMapper>(node, beam);
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, {});
    return channel;
}

template<typename T>
achilles::Channel<achilles::FourVector> BuildChannel(achilles::NuclearModel *model, size_t nlep, size_t nhad,
                                                 std::shared_ptr<achilles::Beam> beam,
                                                 const std::vector<double> &masses) {
    achilles::Channel<achilles::FourVector> channel;
    channel.mapping = achilles::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .FinalState(T::Name(), masses).build();
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    return channel;
}

#ifdef ENABLE_BSM
template<typename T>
achilles::Channel<achilles::FourVector> BuildChannelSherpa(achilles::NuclearModel *model, size_t nlep, size_t nhad,
                                                       std::shared_ptr<achilles::Beam> beam,
                                                       const std::vector<double> &masses) {
    achilles::Channel<achilles::FourVector> channel;
    auto massesGeV = masses;
    for(auto &mass : massesGeV) mass /= (1000*1000);
    channel.mapping = achilles::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .SherpaFinalState(T::Name(), massesGeV).build();
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    return channel;
}

achilles::Channel<achilles::FourVector> BuildGenChannel(achilles::NuclearModel *model, size_t nlep, size_t nhad,
                                                    std::shared_ptr<achilles::Beam> beam,
                                                    std::unique_ptr<PHASIC::Channels> final_state,
                                                    const std::vector<double> &masses) {
    achilles::Channel<achilles::FourVector> channel;
    channel.mapping = achilles::PSBuilder(nlep, nhad).Beam(beam, 1)
                                                   .Hadron(model -> PhaseSpace(), masses)
                                                   .GenFinalState(std::move(final_state)).build();
    achilles::AdaptiveMap map(channel.mapping -> NDims(), 2);
    channel.integrator = achilles::Vegas(map, achilles::VegasParams{});
    return channel;
}
#endif

achilles::EventGen::EventGen(const std::string &configFile,
                             std::vector<std::string> shargs) {
    config = YAML::LoadFile(configFile);

    // Setup random number generator
    auto seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    if(config["Initialize"]["Seed"])
        if(config["Initialize"]["Seed"].as<int>() > 0)
            seed = config["Initialize"]["Seed"].as<unsigned int>();
    spdlog::trace("Seeding generator with: {}", seed);
    Random::Instance().Seed(seed);

    // Load initial state, massess
    spdlog::trace("Initializing the beams");
    beam = std::make_shared<Beam>(config["Beams"].as<Beam>());
    nucleus = std::make_shared<Nucleus>(config["Nucleus"].as<Nucleus>());

    // Set potential for the nucleus
    auto potential_name = config["Nucleus"]["Potential"]["Name"].as<std::string>();
    auto potential = achilles::PotentialFactory::Initialize(potential_name,
                                                          nucleus,
                                                          config["Nucleus"]["Potential"]);
    nucleus -> SetPotential(std::move(potential));
    // Initialize Cascade parameters
    spdlog::debug("Cascade mode: {}", config["Cascade"]["Run"].as<bool>());
    if(config["Cascade"]["Run"].as<bool>()) {
        cascade = std::make_unique<Cascade>(config["Cascade"].as<Cascade>());
    } else {
        cascade = nullptr;
    }

    // Initialize the lepton final states
    spdlog::debug("Initializing the leptonic final states");
    auto leptonicProcess = config["Process"].as<achilles::Process_Info>();
    // TODO: Handle the beam initial state better
    if(beam -> BeamIDs().size() > 1)
        throw std::runtime_error("Multiple processes are not implemented yet. Please use only one beam.");
    leptonicProcess.m_ids.insert(leptonicProcess.m_ids.begin(),
                                 beam -> BeamIDs().begin(), beam -> BeamIDs().end());

    // Initialize the nuclear model
    spdlog::debug("Initializing nuclear model");
    const auto model_name = config["NuclearModel"]["Model"].as<std::string>();
    auto nuclear_model = NuclearModelFactory::Initialize(model_name, config);
    nuclear_model -> AllowedStates(leptonicProcess);
    spdlog::debug("Process: {}", leptonicProcess);

#ifdef ENABLE_BSM
    // Initialize sherpa processes
    achilles::SherpaMEs *sherpa = new achilles::SherpaMEs();
    std::string model = config["Process"]["Model"].as<std::string>();
    std::string param_card = config["Process"]["ParamCard"].as<std::string>();
    if(model == "SM") model = "SM_Nuc";
    shargs.push_back("MODEL=" + model);
    shargs.push_back("UFO_PARAM_CARD=" + param_card);
    sherpa -> Initialize(shargs);
    spdlog::debug("Initializing leptonic currents");
    if(!sherpa -> InitializeProcess(leptonicProcess)) {
        spdlog::error("Cannot initialize hard process");
        exit(1);
    }
    leptonicProcess.m_mom_map = sherpa -> MomentumMap(leptonicProcess.Ids());
#else
    // Dummy call to remove unused error
    shargs.size();
    leptonicProcess.m_mom_map[0] = leptonicProcess.Ids()[0];
    leptonicProcess.m_mom_map[1] = leptonicProcess.Ids()[1];
    leptonicProcess.m_mom_map[2] = leptonicProcess.Ids()[2];
    leptonicProcess.m_mom_map[3] = leptonicProcess.Ids()[3];
#endif

    // Initialize hard cross-sections
    spdlog::debug("Initializing hard interaction");
    scattering = std::make_shared<HardScattering>();
    scattering -> SetProcess(leptonicProcess);
#ifdef ENABLE_BSM
    scattering -> SetSherpa(sherpa);
#endif
    scattering -> SetNuclear(std::move(nuclear_model));

    // Setup channels
    spdlog::debug("Initializing phase space");
    std::vector<double> masses = scattering -> Process().Masses();
    spdlog::trace("Masses = [{}]", fmt::join(masses.begin(), masses.end(), ", "));
    if(config["TestingPS"]) {
        Channel<FourVector> channel = BuildChannelTest(config["TestingPS"], beam);
        integrand.AddChannel(std::move(channel));
    } else {
#ifndef ENABLE_BSM
        if(scattering -> Process().Multiplicity() == 4) {
            Channel<FourVector> channel0 = BuildChannel<TwoBodyMapper>(scattering -> Nuclear(), 2, 2,
                                                                       beam, masses);
            integrand.AddChannel(std::move(channel0));
        } else {
            const std::string error = fmt::format("Leptonic Tensor can only handle 2->2 processes without "
                                                  "BSM being enabled. "
                                                  "Got a 2->{} process", leptonicProcess.m_ids.size());
            throw std::runtime_error(error);
        }
#else
        auto channels = sherpa -> GenerateChannels(scattering -> Process().Ids());
        size_t count = 0;
        for(auto & chan : channels) {
            Channel<FourVector> channel = BuildGenChannel(scattering -> Nuclear(), 
                                                          scattering -> Process().m_ids.size(), 2,
                                                          beam, std::move(chan), masses);
            integrand.AddChannel(std::move(channel));
            spdlog::info("Adding Channel{}", count++);
        }
#endif
    }

    // Setup Multichannel integrator
    // auto params = config["Integration"]["Params"].as<MultiChannelParams>();
    integrator = MultiChannel(integrand.NDims(), integrand.NChannels(), {1000, 2});

    // Decide whether to rotate events to be measured w.r.t. the lepton plane
    if(config["Main"]["DoRotate"])
        doRotate = config["Main"]["DoRotate"].as<bool>();

    // Setup Cuts
    if(config["Main"]["HardCuts"])
        doHardCuts = config["Main"]["HardCuts"].as<bool>();
    spdlog::info("Apply hard cuts? {}", doHardCuts);
    hard_cuts = config["HardCuts"].as<achilles::CutCollection>();

    // if(config["Main"]["EventCuts"])
    //     doEventCuts = config["Main"]["EventCuts"].as<bool>();
    // spdlog::info("Apply event cuts? {}", doEventCuts);
    // event_cuts = config["EventCuts"].as<achilles::CutCollection>();

    // Setup outputs
    auto output = config["Main"]["Output"];
    bool zipped = true;
    if(output["Zipped"])
        zipped = output["Zipped"].as<bool>();
    spdlog::trace("Outputing as {} format", output["Format"].as<std::string>());
    if(output["Format"].as<std::string>() == "Achilles") {
        writer = std::make_unique<AchillesWriter>(output["Name"].as<std::string>(), zipped);
#ifdef ENABLE_HEPMC3
    } else if(output["Format"].as<std::string>() == "HepMC3") {
        writer = std::make_unique<HepMC3Writer>(output["Name"].as<std::string>(), zipped);
#endif
    } else {
        std::string msg = fmt::format("Achilles: Invalid output format requested {}",
                                      output["Format"].as<std::string>());
        throw std::runtime_error(msg);
    }
    writer -> WriteHeader(configFile);
}

void achilles::EventGen::Initialize() {
    // TODO: Clean up loading of previous results
    auto func = [&](const std::vector<FourVector> &mom, const double &wgt) {
        return GenerateEvent(mom, wgt);
    };
    // TODO: Loading the saved data is broken
    // try {
    //     YAML::Node old_results = YAML::LoadFile("results.yml");
    //     integrator = old_results["Multichannel"].as<MultiChannel>();
    //     integrand = old_results["Channels"].as<Integrand<FourVector>>();
    //     YAML::Node results;
    //     results["Multichannel"] = integrator;
    //     results["Channels"] = integrand;
    //     integrand.Function() = func;
    // } catch(const YAML::BadFile &e) {
        spdlog::info("Initializing integrator.");
        integrand.Function() = func;
        if(config["Initialize"]["Accuracy"])
            integrator.Parameters().rtol = config["Initialize"]["Accuracy"].as<double>();
        integrator.Optimize(integrand);
        integrator.Summary();

        YAML::Node results;
        results["Multichannel"] = integrator;
        results["Channels"] = integrand;

        std::ofstream fresults("results.yml");
        fresults << results;
        fresults.close();
    // }
}

double achilles::EventGen::PropogateError(const double p_1, const double p_2, const double error_1, const double error_2) const {
    double num_squared = 2 * pow(p_1, 2) * pow(error_1, 2) + 2 * pow(p_2, 2) * pow(error_2, 2);
    double denom_squared = pow(p_1, 2) + pow(p_2, 2);
    return ((0.5 * sqrt(num_squared)) / sqrt(denom_squared));
}

/* void achilles::EventGen::MakeHist(int arr[], int n){
    int max_elem = *std::max_element(arr, arr + n);
    for (int i = max_elem; i >= 0; i--) {
        std::cout.width(2);
        std::cout << std::right << i << " | ";
        for (int j = 0; j < n; j++) {
            if (arr[j] >= i) {
                std::cout << " x ";
            }
            else {
                std::cout << " ";
            }
        }
        std::cout << "\n";
    }
    for (int i = 0; i < n + 3; i++) {
        std::cout << "---";
        std::cout << "\n";
        std::cout << " ";
    }
    for (int i = 0; i < n; i++) {
        std::cout.width(2);
        std::cout << std::right << arr[i] << " ";
    }
} */

void achilles::EventGen::GenerateEvents() {
    outputEvents = true;
    runCascade = config["Cascade"]["Run"].as<bool>();
    integrator.Parameters().ncalls = config["Main"]["NEvents"].as<size_t>();
    integrator(integrand);
    fmt::print("\n");
    auto result = integrator.Summary();
    fmt::print("Integral = {:^8.5e} +/- {:^8.5e} ({:^8.5e} %)\n",
               result.results.back().Mean(), result.results.back().Error(),
               result.results.back().Error() / result.results.back().Mean()*100);

    /* // set Target_Theta
    std::cout << "Enter target theta (degrees)" << "\n";
    Target_theta = (double)(getchar()); */

    /* // CREATE HISTOGRAM OF Q (MOVED TO GenerateEvent)
    // export data to txt
    // bool anti = ((scattering -> Process().m_ids)[0].AsInt() < 0);
    // if ( ((!anti) && (Amps2[1] == Amps2[1])) || ((anti) && Amps2[0] == Amps2[0]) ) {
    fmt::print("q_0 = {:^8.5e}\n", Q0);
    try {
        std::cout << "Writing data to file" << "\n";
        std::ofstream q_hist("/Users/sherry/Desktop/Fermilab/q_hist_2.txt", std::ios_base::app);
        if (q_hist.is_open()) {
            q_hist << Q0 << "\n";
        }
        else {
            std::cout << "There was a problem opening the file" << "\n";
        }
    }
    catch (const char* msg) {
        std::cerr << msg << "\n";
    }
    std::cout << "Done!\n";*/
    // }
    /* // read file from txt into array
    std::ifstream q_hist_in;
    q_hist_in.open("q_hist.txt");
    double q_data[21];
    double elem;
    if (q_hist_in.is_open())
    {
        int i = 0;
        while (q_hist_in >> elem) {
            q_data[i++] = elem;
        }
    }
    fmt::print("{:^8.5e}\n", q_data[0]);
    fmt::print("{:^8.5e}\n", q_data[10]);
    fmt::print("{:^8.5e}\n", q_data[20]); */


    /* // REPRODUCE FIG 5 PANEL 2 (MOVED TO GenerateEvent)
    bool anti = ((scattering -> Process().m_ids)[0].AsInt() < 0);
    if ((anti) && Amps2[0] == Amps2[0]) {
        // print results
        fmt::print("Beam Energy = {:^8.5e}\n", config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>());
        fmt::print("q_0 = {:^8.5e}\n", Q0);
        fmt::print("Polarization_L (k = 0) = {:^8.5e} +/- {:^8.5e}\n", Polarization_l[0].Mean(), Polarization_l[0].Error());
        // export resuults to file
        try {
            std::cout << "Writing data to file" << "\n";
            std::ofstream fig5_anti_tau_data("/Users/sherry/Desktop/Fermilab/fig5_anti_tau_data.txt", std::ios_base::app);
            if (fig5_anti_tau_data.is_open()) {
                // anti_tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << Polarization_l[0].Mean() << "\t" << Polarization_l[0].Error() << "\n";
                fig5_anti_tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << Q0 << "\t" << Polarization_l[0].Mean() << "\t"
                << Polarization_l[0].Error() << "\n";
            }
            else {
                std::cout << "There was a problem opening the file" << "\n";
            }
        }
        catch (const char* msg) {
            std::cerr << msg << "\n";
        }
        std::cout << "Done!\n";
    }
    else if ((!anti) && (Amps2[1] == Amps2[1]) && (Q0 < 2.1) && (Q0 > 1.9)) {
        // print results
        fmt::print("Beam Energy = {:^8.5e}\n", config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>());
        fmt::print("q_0 = {:^8.5e}\n", Q0);
        fmt::print("Polarization_L (k = 1) = {:^8.5e} +/- {:^8.5e}\n", Polarization_l[1].Mean(), Polarization_l[1].Error());
        // export resuults to file
        try {
            std::cout << "Writing data to file" << "\n";
            std::ofstream fig5_tau_data("/Users/sherry/Desktop/Fermilab/fig5_tau_data.txt", std::ios_base::app);
            if (fig5_tau_data.is_open()) {
                // tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << Polarization_l[1].Mean() << "\t" << Polarization_l[1].Error() << "\n";
                fig5_tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << Q0 << "\t" << Polarization_l[1].Mean() << "\t"
                << Polarization_l[1].Error() << "\n";
            }
            else {
                std::cout << "There was a problem opening the file" << "\n";
            }
        }
        catch (const char* msg) {
            std::cerr << msg << "\n";
        }
        std::cout << "Done!\n";
    } */

    /* // REPRODUCE FIG 7 (WAITING ON JOSH AND NOEMI'S APPROVAL)
    // if PID = 18, anti tau neutrino --> k = 0
    // if PID = 16, tau neutrino --> k = 1
    // TODO: replace with variable
    bool anti = ((scattering -> Process().m_ids)[0].AsInt() < 0);
    // if anti tau neutrino and amps2[0] is not nan
    if ((anti) && Amps2[0] == Amps2[0]) {
        // print results
        fmt::print("Beam Energy = {:^8.5e}\n", config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>());
        fmt::print("Polarization_L (k = 0) = {:^8.5e} +/- {:^8.5e}\n", Polarization_l[0].Mean(), Polarization_l[0].Error());
        fmt::print("Polarization_T (k = 0) = {:^8.5e} +/- {:^8.5e}\n", Polarization_t[0].Mean(), Polarization_t[0].Error());
        fmt::print("Degree of Polarization (k = 0) = {:^8.5e} +/- {:^8.5e}\n", sqrt( pow(Polarization_l[0].Mean(), 2) + pow(Polarization_t[0].Mean(), 2) ),
        PropogateError(Polarization_l[0].Mean(), Polarization_t[0].Mean(), Polarization_l[0].Error(), Polarization_t[0].Error()));
        // export results to file
        try {
            std::cout << "Writing data to file" << "\n";
            std::ofstream anti_tau_data("/Users/sherry/Desktop/Fermilab/fig7_anti_tau_data.txt", std::ios_base::app);
            if (anti_tau_data.is_open()) {
                // anti_tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << Polarization_l[0].Mean() << "\t" << Polarization_l[0].Error() << "\n";
                anti_tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << sqrt( pow(Polarization_l[0].Mean(), 2) + pow(Polarization_t[0].Mean(), 2) )
                << "\t" << PropogateError(Polarization_l[0].Mean(), Polarization_t[0].Mean(), Polarization_l[0].Error(), Polarization_t[0].Error()) << "\n";
            }
            else {
                std::cout << "There was a problem opening the file" << "\n";
            }
        }
        catch (const char* msg) {
            std::cerr << msg << "\n";
        }
        std::cout << "Done!\n";
    }
    // if tau neutrino and Amps2[1] is not nan
    else if ((!anti) && (Amps2[1] == Amps2[1])) {
        // print results
        fmt::print("Beam Energy = {:^8.5e}\n", config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>());
        fmt::print("Polarization_L (k = 1) = {:^8.5e} +/- {:^8.5e}\n", Polarization_l[1].Mean(), Polarization_l[1].Error());
        fmt::print("Polarization_T (k = 1) = {:^8.5e} +/- {:^8.5e}\n", Polarization_t[1].Mean(), Polarization_t[1].Error());
        fmt::print("Degree of Polarization (k = 1) = {:^8.5e} +/- {:^8.5e}\n", sqrt( pow(Polarization_l[1].Mean(), 2) + pow(Polarization_t[1].Mean(), 2) ),
        PropogateError(Polarization_l[1].Mean(), Polarization_t[1].Mean(), Polarization_l[1].Error(), Polarization_t[1].Error()));
        // export results to file
        try {
            std::cout << "Writing data to file" << "\n";
            std::ofstream tau_data("/Users/sherry/Desktop/Fermilab/fig7_tau_data.txt", std::ios_base::app);
            if (tau_data.is_open()) {
                // tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << Polarization_l[1].Mean() << "\t" << Polarization_l[1].Error() << "\n";
                tau_data << config["Beams"][0]["Beam"]["Beam Params"]["Energy"].as<double>() << "\t" << sqrt( pow(Polarization_l[1].Mean(), 2) + pow(Polarization_t[1].Mean(), 2) )
                << "\t" << PropogateError(Polarization_l[1].Mean(), Polarization_t[1].Mean(), Polarization_l[1].Error(), Polarization_t[1].Error()) << "\n";
            }
            else {
                std::cout << "There was a problem opening the file" << "\n";
            }
        }
        catch (const char* msg) {
            std::cerr << msg << "\n";
        }
        std::cout << "Done!\n";
    } */
} 

double achilles::EventGen::GenerateEvent(const std::vector<FourVector> &mom, const double &wgt) {
    if(outputEvents) {
        static size_t ievent = 0;
        constexpr size_t statusUpdate = 10000;
        if(++ievent % statusUpdate == 0) {
            fmt::print("Generated {} / {} events\r", ievent, config["Main"]["NEvents"].as<size_t>());
        }
    }
    // Initialize the event, which generates the nuclear configuration
    // and initializes the beam particle for the event
    Event event(nucleus, mom, wgt);

    // Initialize the particle ids for the processes
    const auto pids = scattering -> Process().m_ids;

    spdlog::debug("Event Phase Space:");
    size_t idx = 0;
    for(const auto &momentum : event.Momentum()) {
        spdlog::debug("\t{}: {} (M2 = {})", ++idx, momentum, momentum.M2());
    }

    // Calculate the hard cross sections and select one for initial state
    spdlog::debug("Calculating cross section");

    // Obtain the cross section for the event 
    auto xsecs = scattering -> CrossSection(event);

    // Initialize the event
    if(!scattering -> FillEvent(event, xsecs)) {
        if(outputEvents) {
            event.SetMEWeight(0);
            spdlog::trace("Outputting the event");
            writer -> Write(event);
        }
        return 0;
    }

    spdlog::trace("Leptons:");
    idx = 0;
    for(const auto &particle : event.Leptons()) {
        spdlog::trace("\t{}: {}", ++idx, particle);
    }

    spdlog::trace("Hadrons:");
    idx = 0;
    for(const auto &particle : event.Hadrons()) {
        spdlog::trace("\t{}: {}", ++idx, particle);
    }

    spdlog::trace("Weight: {}", event.Weight());

    // if((event.Momentum()[3]+event.Momentum()[4]).M() < 400) {
    //     spdlog::info("Mass issue");
    //     spdlog::drop("achilles");
    //     CreateLogger(0, 5);
    // }

    // Perform hard cuts
    if(doHardCuts) {
        spdlog::debug("Making hard cuts");
        if(!MakeCuts(event)) {
            // Short-circuit the evaluation
            // We want Vegas to adapt to avoid these points, i.e.,
            // the integrand should be interpreted as zero in this region
            if(outputEvents) {
                event.SetMEWeight(0);
                writer -> Write(event);
            }
            return 0;
        }
    }

    // Run the cascade if needed
    if(runCascade) {
        spdlog::debug("Runnning cascade");
        cascade -> Evolve(&event);
    } else {
        for(auto & nucleon : event.CurrentNucleus()->Nucleons()) {
            if(nucleon.Status() == ParticleStatus::propagating) {
                nucleon.Status() = ParticleStatus::final_state;
            }
        }
    }

    // Write out events
    if(outputEvents) {
        // Rotate cuts into plane of outgoing electron before writing
        if (doRotate)
            Rotate(event);

        // Perform event-level final cuts before writing
        bool outputCurrentEvent = true;
        // if(doEventCuts){
        //     spdlog::debug("Making event cuts");
        //     outputCurrentEvent = MakeEventCuts(event);
        // }

        if(outputCurrentEvent) {
            // Keep a running total of the number of surviving events
            event.Finalize();
            writer -> Write(event);
            // TODO: generalize to hadronCurrent.size()
            for (size_t k = 0; k < 2; ++k) {
                // update polarization
                Polarization_l[k] += event.get_polarization_l()[k];
                Polarization_t[k] += event.get_polarization_t()[k];

                // update amps2[k]
                Amps2[k] = event.get_amps2()[k];

                // update q_0 and produce q_0 hist data
                Q0 = event.Momentum()[1].E() - event.Momentum().back().E();
                /* fmt::print("q_0 = {:^8.5e}\n", Q0);
                try {
                    std::cout << "Writing data to file" << "\n";
                    std::ofstream q_hist("/Users/sherry/Desktop/Fermilab/q_hist_2.txt", std::ios_base::app);
                    if (q_hist.is_open()) {
                        q_hist << Q0 << "\n";
                    }
                    else {
                        std::cout << "There was a problem opening the file" << "\n";
                    }
                }
                catch (const char* msg) {
                    std::cerr << msg << "\n";
                }
                std::cout << "Done!\n"; */

                // update theta and produce theta hist data
                Theta = event.Momentum().back().Theta();
                /* fmt::print("theta = {:^8.5e}\n", Theta);
                try {
                    std::cout << "Writing data to file" << "\n";
                    std::ofstream theta_hist("/Users/sherry/Desktop/Fermilab/theta_hist_1.txt", std::ios_base::app);
                    if (theta_hist.is_open()) {
                        theta_hist << Theta << "\n";
                    }
                    else {
                        std::cout << "There was a problem opening the file" << "\n";
                    }
                }
                catch (const char* msg) {
                    std::cerr << msg << "\n";
                }
                std::cout << "Done!\n"; */

                // set target theta
                double Target_theta = 16;
                double Theta_degrees = (Theta / M_PI) * 180;

                // filter out unwanted theta values
                if ((Theta_degrees > (Target_theta + 0.5)) || (Theta_degrees < (Target_theta - 0.5))) {

                    // export P_L data for figure 5 (antineutrino)
                    fmt::print("Polarization_L (k = 0) = {:^8.5e}\n", event.get_polarization_l()[0]);
                    try {
                        std::cout << "Writing data to file" << "\n";
                        std::ofstream pl_vs_q_data("/Users/sherry/Desktop/Fermilab/pl_vs_q_2.txt", std::ios_base::app);
                        if (pl_vs_q_data.is_open()) {
                            pl_vs_q_data << Q0 << "\t" << event.get_polarization_l()[0] << "\n";
                        }
                        else {
                            std::cout << "There was a problem opening the file" << "\n";
                        }
                    }
                    catch (const char* msg) {
                        std::cerr << msg << "\n";
                    }
                    std::cout << "Done!\n";

                    // export P_T data for figure 5 (antineutrino)
                    fmt::print("Polarization_T (k = 0) = {:^8.5e}\n", event.get_polarization_t()[0]);
                    try {
                        std::cout << "Writing data to file" << "\n";
                        std::ofstream pt_vs_q_data("/Users/sherry/Desktop/Fermilab/pt_vs_q_2.txt", std::ios_base::app);
                        if (pt_vs_q_data.is_open()) {
                            pt_vs_q_data << Q0 << "\t" << event.get_polarization_t()[0] << "\n";
                        }
                        else {
                            std::cout << "There was a problem opening the file" << "\n";
                        }
                    }
                    catch (const char* msg) {
                        std::cerr << msg << "\n";
                    }
                    std::cout << "Done!\n";
                }
            }
        }
    }

    // Always return the weight when the event passes the initial hard cut.
    // Even if events do not survive the final event-level cuts, Vegas should
    // still interpret the integrand as nonzero in this region.
    return event.Weight();
}

bool achilles::EventGen::MakeCuts(Event &event) {
    return hard_cuts.EvaluateCuts(event.Particles());
}

// TODO: Create Analysis level cuts
/*
bool achilles::EventGen::MakeEventCuts(Event &event) {
    // Run through all particles in the event
    for (const auto& pair : event_cuts) {
        auto pid = pair.first;
        auto cut = pair.second;
        bool pid_passed = false;
        for (const auto& particle : event.Particles()){
            // Restrict to matching final-state particles
            if(particle.IsFinal() && particle.ID() == pid)
                // Keep: at least one particle (of a given PID) survives the cut
                if(cut(particle.Momentum())){
                    pid_passed = true;
                    break;
                }
        }
        // Reject: no particles (of a given PID) satisfy the cut
        if(!pid_passed)
            return false;
    }
    return true;
}*/

void achilles::EventGen::Rotate(Event &event) {
    // Isolate the azimuthal angle of the outgoing electron
    double phi = 0.0;
    for(const auto & particle : event.Particles()){
        if(particle.ID() == PID::electron() && particle.IsFinal()){
            phi = particle.Momentum().Phi();
        }
    }
    // Rotate the coordiantes of particles so that all azimuthal angles phi are
    // measured with respect to the leptonic plane
    std::array<double, 9> rotation = {
        cos(phi),  sin(phi), 0,
        -sin(phi), cos(phi), 0,
        0,         0,        1};
    event.Rotate(rotation);
}
