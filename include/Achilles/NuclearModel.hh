#ifndef NUCLEAR_MODEL_HH
#define NUCLEAR_MODEL_HH

#include "Achilles/Event.hh"
#include "Achilles/Factory.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/SpectralFunction.hh"

#include "yaml-cpp/node/node.h"

#include <complex>

namespace achilles {

class PID;
class PSBuilder;
class Spinor;

enum class NuclearMode {
    None = -1,
    Coherent,
    Quasielastic,
    MesonExchangeCurrent,
    Interference_QE_MEC,
    Resonance,
    ShallowInelastic,
    DeepInelastic
};

class NuclearModel {
    public:
        using Current = std::vector<std::vector<std::complex<double>>>;
        using Currents = std::map<int, Current>;
        using FFInfoMap = std::map<int, std::vector<FormFactorInfo>>;
        using FormFactorMap = std::map<FormFactorInfo::Type, std::complex<double>>;

        NuclearModel() = default;
        NuclearModel(const YAML::Node&, FormFactorBuilder&);
        NuclearModel(const NuclearModel&) = delete;
        NuclearModel(NuclearModel&&) = default;
        NuclearModel& operator=(const NuclearModel&) = delete;
        NuclearModel& operator=(NuclearModel&&) = default;
        virtual ~NuclearModel() = default;

        virtual NuclearMode Mode() const = 0;
        virtual std::string PhaseSpace() const = 0;
        virtual std::vector<Currents> CalcCurrents(const Event&, const std::vector<FFInfoMap>&) const = 0;
        virtual void AllowedStates(Process_Info&) const = 0;
        virtual size_t NSpins() const = 0;
        virtual bool FillNucleus(Event&, const std::vector<double>&) const = 0;

        static std::string Name() { return "Nuclear Model"; }

    protected:
        FormFactor::Values EvalFormFactor(double q2) const { return m_form_factor -> operator()(q2); }
        FormFactorMap CouplingsFF(const FormFactor::Values&,
                                  const std::vector<FormFactorInfo>&) const;
        static YAML::Node LoadFormFactor(const YAML::Node&);

    private:
        std::unique_ptr<FormFactor> m_form_factor{nullptr};
};


template<typename Derived>
using RegistrableNuclearModel = Registrable<NuclearModel, Derived, const YAML::Node&>;
using NuclearModelFactory = Factory<NuclearModel, const YAML::Node&>;

class Coherent : public NuclearModel, RegistrableNuclearModel<Coherent> {
    public:
        Coherent(const YAML::Node&, const YAML::Node&, FormFactorBuilder&);

        NuclearMode Mode() const override { return NuclearMode::Coherent; }
        std::string PhaseSpace() const override { return Name(); }
        std::vector<Currents> CalcCurrents(const Event&, const std::vector<FFInfoMap>&) const override;
        void AllowedStates(Process_Info&) const override;
        size_t NSpins() const override { return 1; }
        bool FillNucleus(Event&, const std::vector<double>&) const override;

        // Required factory methods
        static std::unique_ptr<NuclearModel> Construct(const YAML::Node&);
        static std::string Name() { return "Coherent"; }

    private:
        PID nucleus_pid;
};

class QESpectral : public NuclearModel, RegistrableNuclearModel<QESpectral> {
    public:
        QESpectral(const YAML::Node&, const YAML::Node&, FormFactorBuilder&);

        NuclearMode Mode() const override { return NuclearMode::Quasielastic; }
        std::string PhaseSpace() const override { return Name(); }
        std::vector<Currents> CalcCurrents(const Event&, const std::vector<FFInfoMap>&) const override;
        void AllowedStates(Process_Info&) const override;
        size_t NSpins() const override { return 4; }
        bool FillNucleus(Event&, const std::vector<double>&) const override;

        // Required factory methods
        static std::unique_ptr<NuclearModel> Construct(const YAML::Node&);
        static std::string Name() { return "QESpectral"; }

    private:
        bool b_ward{};
        Current HadronicCurrent(const std::array<Spinor, 2>&, const std::array<Spinor, 2>&,
                                const FourVector&, const FormFactorMap&) const;
        SpectralFunction spectral_proton, spectral_neutron; 
};

}

#endif
