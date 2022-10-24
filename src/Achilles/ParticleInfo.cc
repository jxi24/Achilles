#include "Achilles/ParticleInfo.hh"
#include "Achilles/Exception.hh"

#include <iostream>
#include <memory>

using achilles::PID;
using achilles::ParticleInfoEntry;
using achilles::ParticleInfo;

// achilles::Database::ParticleDB achilles::Database::particleDB{};

// PID::PID(const long int &id_, bool validate) : id(PID::undefined()) {
//     if(!validate)
//          id = std::abs(id_);
//     else
//          valid(std::abs(id_));
//
//     id = std::abs(id_);
//     // Check if fundamental SM particles
//     if(id > 0 && id < 30) {
//     }
// }

void achilles::ParticleInfo::BuildDatabase(const std::string &datafile) {
    YAML::Node particleYAML = YAML::LoadFile(datafile);
    auto particles = particleYAML["Particles"];
    for(auto particle : particles) {
        auto entry = std::make_shared<ParticleInfoEntry>(particle["Particle"].as<ParticleInfoEntry>());
        particleDB.emplace(entry -> id, entry);
        nameToPID.emplace(entry -> idname, entry -> id);
    }
    PrintDatabase();
}

void achilles::ParticleInfo::PrintDatabase() {
    fmt::print("{:>10s} {:<20s} {:<20s} {:^10s}    {:^10s}\n",
               "PID", "Name", "Anti-name", "Mass (MeV)", "Width (MeV)");
    for(const auto &part : particleDB) {
        std::cout <<  *(part.second) << "\n";
    }
}

ParticleInfoEntry::ParticleInfoEntry(const achilles::PID &id_, const double &mass_,
                                     const double &width_, const int &icharge_,
                                     const int &strong_, const int &spin_,
                                     const int &stable_, const int &majorana_, 
                                     const bool &massive_, const bool &hadron_,
                                     std::string idname_, std::string antiname_) 
    : id(id_), mass(mass_), hmass(mass_), width(width_), icharge(icharge_), strong(strong_),
      spin(spin_), stable(stable_), majorana(majorana_), massive(massive_), hadron(hadron_), 
      idname(std::move(idname_)), antiname(std::move(antiname_)) {}

std::ostream& operator<<(std::ostream &os, const ParticleInfoEntry &entry) {
    os << entry.ToString();
    return os;
}

ParticleInfo::ParticleDB ParticleInfo::particleDB;
std::map<std::string, achilles::PID> ParticleInfo::nameToPID;

bool ParticleInfo::IsBaryon() const noexcept {
    if(IntID() % 10000 < 1000) return false;
    return true;
}

bool ParticleInfo::IsBHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID()-100*IntID()/100 < 10) return false;
    if((IntID()-100*IntID()/100) / 10 == 5) return true;
    if((IntID()-1000*IntID()/1000) / 100 == 5) return true;
    if((IntID()-10000*IntID()/10000) / 1000 == 5) return true;
    return false;
}

bool ParticleInfo::IsCHadron() const noexcept {
    if(IntID() < 100) return false;
    if(IntID() - 100*IntID()/100 < 10) return false;
    if((IntID()-100*IntID()/100) / 10 == 4) return true;
    if((IntID()-1000*IntID()/1000) / 100 == 4) return true;
    if((IntID()-10000*IntID()/10000) / 1000 == 4) return true;
    return false;
}

double ParticleInfo::GenerateLifeTime() const {
    throw not_implemented_error("Lifetimes for particles are not implemented");
    return 0.0;
}

bool ParticleInfo::IsStable() const noexcept {
    if(info -> stable == 0) return false;
    if(info -> stable == 1) return true;
    if(info -> stable == 2 && !IsAnti()) return true;
    if(info -> stable == 3 && IsAnti()) return true;
    return false;
}
