#include <cmath>
#include <map>
#include <regex>
#include <fstream>
#include <iostream>

#include <cmath>
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "Achilles/Constants.hh"
#include "Achilles/Exception.hh"
#include "Achilles/ThreeVector.hh"
#include "Achilles/Particle.hh"
#include "Achilles/Nucleus.hh"
#include "Achilles/Utilities.hh"

using namespace achilles;

const std::map<std::size_t, std::string> Nucleus::ZToName = {
    {0, "mfp"},
    {1, "H"},
    {2, "He"},
    {3, "Li"},
    {6, "C"},
    {8, "O"},
    {13, "Al"},
    {18, "Ar"},
    {20, "Ca"},
    {26, "Fe"},
};

Nucleus::Nucleus(const std::size_t& Z, const std::size_t& A, const double& bEnergy,
                 const double& kf, const std::string& densityFilename, const FermiGasType& fgType,
                 std::unique_ptr<Density> _density) 
                        : binding(bEnergy), fermiMomentum(kf), fermiGas(fgType),
                          density(std::move(_density)) {
    
    if(Z > A) {
        std::string errorMsg = "Requires the number of protons to be less than the total";
        errorMsg += " number of nucleons. Got " + std::to_string(Z);
        errorMsg += " protons and " + std::to_string(A) + " nucleons";
        throw nucleus_error(errorMsg);
    }
    
    nucleons.resize(A);
    protons.resize(Z);
    neutrons.resize(A-Z);
    // TODO: Refactor elsewhere in the code, maybe make dynamic?
    // spdlog::info("Nucleus: inferring nuclear radius using 0.16 nucleons/fm^3.");
    // constexpr double nucDensity = 0.16;
    // radius = std::cbrt(static_cast<double>(A) / (4 / 3 * M_PI * nucDensity));

    std::ifstream densityFile(densityFilename);
    if(!densityFile.is_open())
        throw nucleus_error(fmt::format("Density file {} does not exist.", densityFilename));
    std::string lineContent;
   
    constexpr size_t HeaderLength = 16;
    for(size_t i = 0; i < HeaderLength; ++i) {       
        std::getline(densityFile, lineContent);
    }

    double radius_{}, density_{}, densityErr{};
    std::vector<double> vecRadius, vecDensity;
    constexpr double minDensity = 1E-6;
    while(densityFile >> radius_ >> density_ >> densityErr) {
        if(density_ < minDensity && radius == 0) radius=radius_;
        vecRadius.push_back(std::move(radius_));
        vecDensity.push_back(std::move(density_));
    }

    rhoInterp.SetData(vecRadius, vecDensity);
    rhoInterp.CubicSpline();
    
    // Ensure the number of protons and neutrons are correct
    // NOTE: This only is checked at startup, so if density returns a varying number of nucleons it will 
    // not necessarily be caught 
    auto particles = density -> GetConfiguration();
    if(particles.size() != nucleons.size())
        throw nucleus_error("Invalid density function! Incorrect number of nucleons.");

    std::size_t nProtons = 0, nNeutrons = 0;
    for(auto particle : particles) {
        if(particle.ID() == PID::proton()) nProtons++;
        if(particle.ID() == PID::neutron()) nNeutrons++;
    }

    if(nProtons != NProtons() || nNeutrons != NNeutrons())
        throw nucleus_error("Invalid density function! Incorrect number of protons or neutrons.");

    static constexpr int IDBase = 1000000000;
    static constexpr int ZBase = 10000;
    static constexpr int ABase = 10;
    int ID = IDBase + ZBase*static_cast<int>(NProtons()) + ABase*static_cast<int>(NNucleons());
    m_pid = PID{ID};
}

// achilles::PID Nucleus::ID() const {
//     // Output format based on PDG Monte-Carlo PIDs
//     // Nuclear codes are given as a 10 digit number:
//     // +/- 10LZZZAAAI
//     // L: number of strange baryons
//     // Z: number of protons
//     // A: number of nucleons
//     // I: excited state (0 is ground state)
//     static constexpr int IDBase = 1000000000;
//     static constexpr int ZBase = 10000;
//     static constexpr int ABase = 10;
//     int ID = IDBase + ZBase*static_cast<int>(NProtons()) + ABase*static_cast<int>(NNucleons());
//     return PID{ID};
// }

void Nucleus::SetNucleons(Particles& _nucleons) noexcept {
    nucleons = _nucleons;
    std::size_t idx = 0;
    std::size_t proton_idx = 0;
    std::size_t neutron_idx = 0;
    for(auto particle : nucleons) {
        if(particle.ID() == PID::proton()) {
            if(proton_idx >= protons.size()) {
                protons.push_back(particle);
                proton_idx++;
            } else protons[proton_idx++] = particle;
            protonLoc.push_back(idx++);
        }
        else if(particle.ID() == PID::neutron()) {
            if(neutron_idx >= neutrons.size()) {
                neutrons.push_back(particle);
                neutron_idx++;
            } else neutrons[neutron_idx++] = particle;
            neutronLoc.push_back(idx++);
        }
    }
}

void Nucleus::GenerateConfig() {
    // Get a configuration from the density function
    Particles particles = density -> GetConfiguration();

    for(Particle& particle : particles) {
        // Set momentum for each nucleon
        auto mom3 = GenerateMomentum(particle.Position().Magnitude());
        double energy2 = pow(particle.Info().Mass(), 2); // Constant::mN*Constant::mN;
        for(auto mom : mom3) energy2 += mom*mom;
        particle.Momentum() = FourVector(sqrt(energy2), mom3[0], mom3[1], mom3[2]);

        // Ensure status is set to background
        particle.Status() = ParticleStatus::background;
    }

    // Update the nucleons in the nucleus
    SetNucleons(particles);
}

const std::array<double, 3> Nucleus::GenerateMomentum(const double &position) noexcept {
    std::array<double, 3> momentum{};
    momentum[0] = Random::Instance().Uniform(0.0,FermiMomentum(position));
    momentum[1] = std::acos(Random::Instance().Uniform(-1.0, 1.0));
    momentum[2] = Random::Instance().Uniform(0.0, 2*M_PI);

    return ToCartesian(momentum);
}

Nucleus Nucleus::MakeNucleus(const std::string& name, const double& bEnergy,
                             const double& fermiMomentum,
                             const std::string& densityFilename, const FermiGasType& fg_type,
                             std::unique_ptr<Density> density) {
    const std::regex regex("([0-9]+)([a-zA-Z]+)");
    std::smatch match;

    if(std::regex_match(name, match, regex)) {
        const std::size_t nucleons = std::stoul(match[1].str());
        const std::size_t protons = NameToZ(match[2].str()); 
        spdlog::info(
            "Nucleus: parsing nuclear name '{0}', expecting a density "
            "with A={1} total nucleons and Z={2} protons.", 
            name, nucleons, protons);
        return Nucleus(protons, nucleons, bEnergy, fermiMomentum, densityFilename,
                       fg_type, std::move(density));
    }

    throw nucleus_error(fmt::format("Invalid nucleus {}.", name));
}

std::size_t Nucleus::NameToZ(const std::string& name) {
    auto it = std::find_if(ZToName.begin(), ZToName.end(),
                           [&name](const std::pair<int, std::string> &p) {
                               return p.second == name;
                          });
    if(it == ZToName.end()) 
        throw nucleus_error(fmt::format("Invalid nucleus: {} does not exist.", name));
    return it -> first;
}

const std::string Nucleus::ToString() const noexcept {
    return std::to_string(NNucleons()) + ZToName.at(NProtons());
}

double Nucleus::FermiMomentum(const double &position) const noexcept { 
    double rho = Rho(position);
    double result{};
    switch(fermiGas) {
        case FermiGasType::Local:
            result = std::cbrt(rho*3*M_PI*M_PI)*Constant::HBARC;
            break;
        case FermiGasType::Global:
            static constexpr double small = 1E-2;
            result = rho < small ? small : fermiMomentum;
            break;
    }

    return result;
}
