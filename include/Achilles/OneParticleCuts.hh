#ifndef ONE_PARTICLE_CUT_HH
#define ONE_PARTICLE_CUT_HH

#include "Achilles/Cuts.hh"

namespace achilles {

class OneParticleCut : public CutBase<double> {
    public:
        // Calculate cut for input FourVector
        virtual bool MakeCut(const FourVector&) const = 0;
        OneParticleCut(const YAML::Node &node) : CutBase(node) {}
        virtual ~OneParticleCut() = default;
        static std::string Name() { return "Single Particle"; }
};

#define ONE_PARTICLE_CUT(CutName) \
    class CutName##Cut : public OneParticleCut, RegistrableCut<OneParticleCut, CutName##Cut> { \
        public: \
            CutName##Cut(const YAML::Node &node) : OneParticleCut(node) {} \
            static std::string Name() { return #CutName; } \
            static std::unique_ptr<OneParticleCut> Construct(const YAML::Node &node) { \
                return std::make_unique<CutName##Cut>(node); \
            } \
            bool MakeCut(const FourVector&) const override;\
    }

ONE_PARTICLE_CUT(Energy);
ONE_PARTICLE_CUT(Momentum);
ONE_PARTICLE_CUT(AngleTheta);
ONE_PARTICLE_CUT(TransverseMomentum);
ONE_PARTICLE_CUT(ETheta2);

class Q0_QECut : public OneParticleCut, RegistrableCut<OneParticleCut, Q0_QECut> {
    public:
        Q0_QECut(const YAML::Node &node) : OneParticleCut(node) {
            eb = node["Eb"].as<double>();
        }
        bool MakeCut(const FourVector&) const;
        static std::string Name() { return "Q0_QE"; }
        static std::unique_ptr<OneParticleCut> Construct(const YAML::Node &node) {
            return std::make_unique<Q0_QECut>(node);
        }

    private:
        double eb;
};

// class EnergyCut : public OneParticleCut, RegistrableCut<OneParticleCut, EnergyCut> {
//     public:
//         EnergyCut(const YAML::Node &node) : OneParticleCut(node) {}
// 
//         static std::string Name() { return "Energy"; } 
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node &node) {
//             return std::make_unique<EnergyCut>(node);
//         }
// 
//         bool MakeCut(const FourVector&) const override { return true; }
// 
//     private:
//         cut_range m_range;
// };

// class MomentumCut : public OneParticleCut, RegistrableCut<OneParticleCut, MomentumCut> {
//     public:
//         MomentumCut() = default;
// 
//         static std::string Name() { return "Momentum"; } 
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<MomentumCut>();
//         }
// 
//         bool MakeCut(const FourVector&) const override { return true; }
// 
//     private:
//         cut_range m_range;
// };
// 
// class AngleThetaCut : public OneParticleCut, RegistrableCut<OneParticleCut, AngleThetaCut> {
//     public:
//         AngleThetaCut() = default;
// 
//         static std::string Name() { return "AngleTheta"; } 
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<AngleThetaCut>();
//         }
// 
//         bool MakeCut(const FourVector&) const override { return true; }
// 
//     private:
//         cut_range m_range;
// };
// 
// class TransverseMomentumCut : public OneParticleCut, RegistrableCut<OneParticleCut, TransverseMomentumCut> {
//     public:
//         TransverseMomentumCut() = default;
// 
//         static std::string Name() { return "TransverseMomentum"; } 
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<TransverseMomentumCut>();
//         }
// 
//         bool MakeCut(const FourVector&) const override { return true; }
// 
//     private:
//         cut_range m_range;
// };
// 
// class ETheta2Cut : public OneParticleCut, RegistrableCut<OneParticleCut, ETheta2Cut> {
//     public:
//         ETheta2Cut() = default;
// 
//         static std::string Name() { return "ETheta2"; } 
//         static std::unique_ptr<OneParticleCut> Construct(const YAML::Node&) {
//             return std::make_unique<ETheta2Cut>();
//         }
// 
//         bool MakeCut(const FourVector&) const override { return true; }
// 
//     private:
//         cut_range m_range;
// };

}

#endif
