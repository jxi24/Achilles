#ifndef CUTS_HH
#define CUTS_HH

#include <utility>

#include "Achilles/ParticleInfo.hh"
#include "fmt/ranges.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "yaml-cpp/yaml.h"
#pragma GCC diagnostic pop

namespace achilles {

class FourVector;

template<class T>
class CutBase {
    public:
        using cut_range = std::pair<T, T>;
        using cut_ranges = std::vector<std::pair<T, T>>;
        CutBase() = default;

        CutBase(const YAML::Node &node) {
            if(node["min"] && node["max"] && node["range"]) {
                throw std::logic_error("CutRange: Invalid syntax. Got min, max, and range");
            } else if(node["min"] || node["max"]) {
                T min = std::numeric_limits<T>::lowest();
                T max = std::numeric_limits<T>::max();
                if(node["min"]) {
                    if(!node["min"].IsScalar())
                        throw std::logic_error("CutRange: Invalid min value");
                    min = node["min"].as<T>();
                }
                if(node["max"]) {
                    if(!node["max"].IsScalar())
                        throw std::logic_error("CutRange: Invalid max value");
                    max = node["max"].as<T>();
                }
                m_range = {{min, max}}; 
            } else if(node["range"]) {
                m_range = node["range"].as<cut_ranges>();
            } else {
                throw std::logic_error("CutRange: Invalid syntax. Missing cut values");
            }
            if(m_range.size() == 1) {
                spdlog::trace("Found cut range: [{}, {}]", m_range[0].first, m_range[0].second);
            } else {
                std::string ranges{};
                for(const auto &range : m_range) {
                    ranges += fmt::format("[{}, {}], ", range.first, range.second);
                }
                spdlog::trace("Found cut range: [{}]", ranges.substr(0, ranges.size()-2));
            }
        }

    protected:
        bool CheckCut(const T &val) const {
            bool result = false;
            // Run through possible cut ranges, keeping the event if a single cut range
            // is satisfied. Note that these cut ranges are all for a single variable,
            // e.g., theta in [0,1] OR [2,3] OR ...
            // Return true if any single cut range is satisfied.
            for(const auto &cut : m_range) {
                result |= (val > cut.first && val < cut.second);
            }

            return result;
        }

    private:
        cut_ranges m_range;
};

template<class Base>
class CutFactory {
    using Constructor = std::function<std::unique_ptr<Base>(const YAML::Node&)>;

    static std::map<std::string, Constructor>& Registry() {
        static std::map<std::string, Constructor> registry;
        return registry;
    }

    public:
        static std::unique_ptr<Base> InitializeCut(const std::string &name, const YAML::Node &node) {
            auto constructor = Registry().at(name);
            return constructor(node);
        }

        template<class Derived>
        static void Register(std::string name) {
            if(IsRegistered(name))
                spdlog::error("{} is already registered!", name);
            spdlog::trace("Registering {} Cut", name);
            Registry()[name] = Derived::Construct;
        }

        static bool IsRegistered(std::string name) {
            return Registry().find(name) != Registry().end();
        }

        static void DisplayCuts() {
            fmt::print("Registered {} cuts:\n", Base::Name());
            for(const auto &registered : Registry()) 
                fmt::print("  - {}\n", registered.first);
        }
};

template<class Base, class Derived>
class RegistrableCut {
    protected:
        RegistrableCut() = default;
        virtual ~RegistrableCut()  {
            if(!m_registered)
                spdlog::error("Error registering cut");
        }

        static bool Register() {
            CutFactory<Base>::template Register<Derived>(Derived::Name());
            return true;
        }

    private:
        static const bool m_registered;
};
template<class Base, class Derived>
const bool RegistrableCut<Base, Derived>::m_registered = RegistrableCut<Base, Derived>::Register();

}

namespace YAML {

template<typename T>
struct convert<std::pair<T, T>> {
    static bool decode(const Node &node, std::pair<T, T> &range) {
        if(node[0].IsScalar() && node[1].IsScalar() && node.size() == 2) {
            range = {node[0].as<T>(), node[1].as<T>()};
            return true;
        }
        return false;
    }
};

template<typename T>
struct convert<std::vector<std::pair<T, T>>> {
    static bool decode(const Node &node, std::vector<std::pair<T, T>> &cutRange) {
        if(!node.IsSequence())
            return false;
        if(node[0].IsScalar() && node[1].IsScalar() && node.size() == 2) {
            cutRange = {{node[0].as<double>(), node[1].as<T>()}};
        } else {
            for(const auto &subNode : node) {
                cutRange.push_back(subNode.as<std::pair<T, T>>());
            }
        }
        return true;
    }
};

}

#endif
