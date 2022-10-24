#ifndef PARTICLE_ID_HH
#define PARTICLE_ID_HH

// The classes in this file are inspired from the implementation found in Sherpa

#include <map>
#include <memory>
#include <spdlog/spdlog.h>
#include <string>
#include <utility>
#include <functional>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {
    class ParticleInfoEntry;
    // class ParticleInfo;
}

std::ostream& operator<<(std::ostream&, const achilles::ParticleInfoEntry&);
// std::ostream& operator<<(std::ostream&, const achilles::ParticleInfo&);

namespace achilles {
    class PID {
        public:
            constexpr PID(const long int &id_ = 0) : id(id_) {}
            // PID(const long int&, bool);
            constexpr long int AsInt() const { return id; }
            constexpr PID Anti() const { return PID{ -id }; }
            constexpr bool operator==(const PID &other) const { return id == other.id; }
            constexpr bool operator!=(const PID &other) const { return id != other.id; }
            constexpr bool operator<(const PID &other) const { return id < other.id; }
            constexpr bool operator>(const PID &other) const { return id > other.id; }
            constexpr PID operator-() const { return PID{ -id }; }
            PID Abs() const { return std::abs(id); }
            constexpr operator long int() const { return id; }
            constexpr operator int() const { return static_cast<int>(id); }

            // Names for common particles
            // undefined
            static constexpr PID undefined() { return PID{ 0 }; }
            // Quarks
            static constexpr PID down() { return PID{ 1 }; } 
            static constexpr PID up() { return PID{ 2 }; } 
            static constexpr PID strange() { return PID{ 3 }; }
            static constexpr PID charm() { return PID{ 4 }; }
            static constexpr PID bottom() { return PID{ 5 }; }
            static constexpr PID top() { return PID{ 6 }; }
            // Leptons
            static constexpr PID electron() { return PID{ 11 }; }
            static constexpr PID nu_electron() { return PID{ 12 }; }
            static constexpr PID muon() { return PID{ 13 }; }
            static constexpr PID nu_muon() { return PID{ 14 }; }
            static constexpr PID tau() { return PID{ 15 }; }
            static constexpr PID nu_tau() { return PID{ 16 }; }
            // Gauge Bosons
            static constexpr PID gluon() { return PID{ 21 }; }
            static constexpr PID photon() { return PID{ 22 }; }
            static constexpr PID Zboson() { return PID{ 23 }; }
            static constexpr PID Wboson() { return PID{ 24 }; }
            static constexpr PID Higgs() { return PID{ 25 }; }
            // Mesons
            static constexpr PID pion0() { return PID{ 111 }; }
            static constexpr PID pionp() { return PID{ 211 }; }
            // Baryons
            static constexpr PID proton() { return PID{ 2212 }; }
            static constexpr PID neutron() { return PID{ 2112 }; }
            // Dummy hadron
            static constexpr PID dummyHadron() { return PID{ 2212 }; }
            // Common Elements
            static constexpr PID dummyNucleus() { return PID{ 1000000000 }; }
            static constexpr PID carbon() { return PID{ 1000060120 }; }
            static constexpr PID argon() { return PID{ 1000180400 }; }

            // Stream Operator
            template<typename OStream>
            friend OStream& operator<<(OStream &os, const PID &pid) {
                os << pid.id;
                return os;
            }

        private:
            // Ensure id matches the numbering scheme defined at:
            // http://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
            // TODO: Implement this function, may not be explicitly needed though
            // bool valid(const long int&);
            long int id;
    };

    class ParticleInfo;

    class ParticleInfoEntry {
        public:
            ParticleInfoEntry()
                : id(PID::undefined()) {}
            ParticleInfoEntry(const PID&, const double&, const double&, const int&, const int&,
                              const int&, const int&, const int&, const bool&, const bool&,
                              std::string, std::string);

            ParticleInfoEntry(const ParticleInfoEntry &info) = default;
            ParticleInfoEntry(ParticleInfoEntry &&info) = default;

            ~ParticleInfoEntry() = default;
            ParticleInfoEntry& operator=(const ParticleInfoEntry &other) = default;
            ParticleInfoEntry& operator=(ParticleInfoEntry &&other) = default;

            bool operator==(const ParticleInfoEntry &other) const noexcept{
                return id == other.id && mass == other.mass && hmass == other.hmass && width == other.width
                    && icharge == other.icharge && strong == other.strong && spin == other.spin && stable == other.stable
                    && majorana == other.majorana && massive == other.massive && hadron == other.hadron && idname == other.idname
                    && antiname == other.antiname;
            }
            bool operator!=(const ParticleInfoEntry &other) const noexcept { 
                return !(*this == other);
            }

            std::string ToString() const noexcept {
                return fmt::format("{:10d} {:20s} {:20s} {: .4e}    {:.4e}", 
                                static_cast<int>(id), idname,
                                antiname, mass, width);
            }

            friend class ParticleInfo;
            friend struct YAML::convert<ParticleInfoEntry>;
            friend std::ostream& ::operator<<(std::ostream&, const ParticleInfoEntry&);

        private:
            PID id;
            double mass{}, hmass{}, width{};
            int icharge{}, strong{};
            int spin{}, stable{}, majorana{};
            bool massive{};
            bool hadron{};
            std::string idname, antiname;
    };

    class ParticleInfo {
        private:
            using ParticleDB = std::map<PID, std::shared_ptr<ParticleInfoEntry>>;
            static ParticleDB particleDB;
            static std::map<std::string, PID> nameToPID; 
            static void BuildDatabase(const std::string&);

        public:
            ParticleInfo(std::shared_ptr<ParticleInfoEntry> info_, const bool &anti_=false)
                : info(std::move(info_)), anti(anti_) {
                InitDatabase("data/Particles.yml");
                if(particleDB.find(info -> id) == particleDB.end())
                    particleDB[info -> id] = info;
                if(anti && info -> majorana == 0) anti = anti_;
            }

            explicit ParticleInfo(const long int &id) : info(nullptr), anti(false) {
                InitDatabase("data/Particles.yml");
                auto it(particleDB.find(static_cast<PID>(std::abs(id))));
                if(it != particleDB.end()) 
                    info = it -> second;
                else
                    throw std::logic_error(fmt::format("Invalid PID: id={}", id));
                if(id < 0 && info -> majorana == 0) anti = true;
            }

            ParticleInfo(PID id, const bool &anti_=false) : info(nullptr), anti(anti_) {
                InitDatabase("data/Particles.yml");
                if(id < PID::undefined()) {
                    id = -id;
                    anti = true;
                }
                auto it(particleDB.find(id));
                if(it == particleDB.end())
                    throw std::logic_error(fmt::format("Invalid PID: id={}", int(id)));
                info = it -> second;
                if(anti_ && info -> majorana == 0) anti = anti_;
            }

            ParticleInfo(const ParticleInfo&) = default;
            ParticleInfo(ParticleInfo&&) = default;
            ~ParticleInfo() = default;
            ParticleInfo& operator=(const ParticleInfo &other) = default;
            ParticleInfo& operator=(ParticleInfo &&other) = default;

            ParticleInfo Anti() { return ParticleInfo(-IntID()); }

            // Property functions
            std::string Name() const noexcept { return anti ? info -> antiname : info -> idname; }
            PID ID() const noexcept { return info -> id; }
            int IntID() const noexcept { return anti ? -static_cast<int>(info->id) : static_cast<int>(info->id); }
            bool IsBaryon() const noexcept;
            bool IsHadron() const noexcept { return info -> hadron; }
            bool IsBHadron() const noexcept;
            bool IsCHadron() const noexcept;
            bool IsAnti() const noexcept { return anti; }
            bool IsFermion() const noexcept { return IntSpin()%2 == 1; }
            bool IsBoson() const noexcept { return IntSpin()%2 == 0; }
            bool IsScalar() const noexcept { return IntSpin() == 0; }
            bool IsVector() const noexcept { return IntSpin() == 2; }
            bool IsTensor() const noexcept { return IntSpin() == 4; }
            bool IsPhoton() const noexcept { return info -> id == PID::photon(); }
            bool IsLepton() const noexcept { return std::abs(IntID()) > 10 && std::abs(IntID()) < 19; }
            bool IsQuark() const noexcept { return IntID() < 10; }
            bool IsGluon() const noexcept { return info -> id == PID::gluon(); }
            bool IsNeutrino() const noexcept { return std::abs(IntID()) == 12 
                                                   || std::abs(IntID()) == 14
                                                   || std::abs(IntID()) == 16; }
            bool IsNucleus() const noexcept { return std::abs(IntID()) > 1000000000; }

            int IntCharge() const noexcept { 
                int charge(info -> icharge); 
                return anti ? -charge : charge;
            }
            double Charge() const noexcept { return static_cast<double>(IntCharge()) / 3; }
            int IntSpin() const noexcept { return info -> spin; }
            double Spin() const noexcept { return static_cast<double>(info -> spin) / 2; }
            bool SelfAnti() const noexcept { return info -> majorana != 0; } 
            bool Majorana() const noexcept { return info -> majorana == 1; }
            int Stable() const noexcept { return info -> stable; }
            bool IsStable() const noexcept;
            bool IsMassive() const noexcept { return info -> mass != 0 ? info -> massive : false; } 
            double Mass() const noexcept { return info -> massive ? info -> mass : 0.0; }
            double Width() const noexcept { return info -> width; }

            double GenerateLifeTime() const;

            bool operator==(const ParticleInfo &other) const noexcept {
                return info == other.info && anti == other.anti;
            }
            bool operator!=(const ParticleInfo &other) const noexcept { return !(*this == other); }

            static ParticleDB& Database() { return particleDB; }
            static void InitDatabase(const std::string &filename) {
                if(particleDB.size() == 0) {
                    particleDB[PID::undefined()] =  std::make_shared<ParticleInfoEntry>(ParticleInfoEntry());
                    BuildDatabase(filename);
                }
            }
            static void PrintDatabase();
            static const std::map<std::string, PID>& NameToPID() { return nameToPID; }

        private:
            std::shared_ptr<ParticleInfoEntry> info;
            bool anti;
    };

    namespace Database {
        static constexpr auto Particle = ParticleInfo::Database;
        static constexpr auto InitParticle = ParticleInfo::InitDatabase;
        static constexpr auto PrintParticle = ParticleInfo::PrintDatabase;
    }

}

namespace YAML {

template<>
struct convert<achilles::PID> {
    static bool decode(const Node &node, achilles::PID &pid) {
        if(node.IsScalar()) {
            pid = achilles::PID(node.as<int>());
            return true;
        }
        return false;
    }
};

template<>
struct convert<std::set<achilles::PID>> {
    static bool decode(const Node &node, std::set<achilles::PID> &pids) {
        if(node.IsScalar()) {
            pids.insert(node.as<achilles::PID>());
            return true;
        }

        for(const auto &subnode : node) {
            pids.insert(subnode.as<achilles::PID>());
        }
        return true;
    }
};

template<>
struct convert<achilles::ParticleInfoEntry> {
    static Node encode(const achilles::ParticleInfoEntry &partInfo) {
        Node node;
        node.push_back(static_cast<int>(partInfo.id));
        node.push_back(partInfo.mass);
        node.push_back(partInfo.width);
        node.push_back(partInfo.icharge);
        node.push_back(partInfo.strong);
        node.push_back(partInfo.spin);
        node.push_back(partInfo.stable);
        node.push_back(partInfo.majorana);
        node.push_back(partInfo.massive);
        node.push_back(partInfo.hadron);
        node.push_back(partInfo.idname);
        node.push_back(partInfo.antiname);

        return node;
    }

    static bool decode(const Node &node, achilles::ParticleInfoEntry &partInfo) {
        if(!node.IsSequence() || node.size() != 12) {
            return false;
        } 

        partInfo.id = static_cast<achilles::PID>(node[0].as<int>());
        partInfo.mass = node[1].as<double>();
        partInfo.width = node[2].as<double>();
        partInfo.icharge = node[3].as<int>();
        partInfo.strong = node[4].as<int>();
        partInfo.spin = node[5].as<int>();
        partInfo.stable = node[6].as<int>();
        partInfo.majorana = node[7].as<int>();
        partInfo.massive = node[8].as<bool>();
        partInfo.hadron = node[9].as<bool>();
        partInfo.idname = node[10].as<std::string>();
        partInfo.antiname = node[11].as<std::string>();

        return true;
    }
};

}

namespace fmt {

template<>
struct formatter<achilles::PID> {
    template<typename ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const achilles::PID &pid, FormatContext &ctx) {
        return format_to(ctx.out(), "{}", pid.AsInt());
    }
};

}

#endif
