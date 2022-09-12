#include "plugins/Sherpa/SherpaMEs.hh"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "COMIX/Main/Single_Process.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Color_Calculator.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "Achilles/Logging.hh"
#include "Achilles/ProcessInfo.hh"
#include "Achilles/ParticleInfo.hh"
#include "Achilles/FormFactor.hh"
#include "Achilles/Utilities.hh"
#include "plugins/Sherpa/Channels.hh"
#include <bitset>

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace achilles;

PHASIC::NLOTypeStringProcessMap_Map m_pmap;
PHASIC::Process_Vector m_procs;

SherpaMEs::~SherpaMEs()
{
  if (p_sherpa) delete p_sherpa;
  for (size_t i(0);i<m_procs.size();++i) delete m_procs[i];
  if (m_pmap.find(nlo_type::lo)!=m_pmap.end()) delete m_pmap[nlo_type::lo];
}

void SherpaMEs::addParameter
(std::vector<char*>& argv,const std::string& val) const
{
  spdlog::info("ShowerMEsSherpa::init: Setting parameter '{}'",val);
  argv.push_back(new char[val.length()+1]);
  strcpy(argv.back(),val.c_str());
}

int SherpaMEs::SherpaVerbosity(int loglevel) const
{
  switch(loglevel) {
  case 0/*trace*/: return 47;
  case 1/*debug*/: return 15;
  case 2/*info*/: return 2;
  case 3/*warn*/: return 1;
  case 4/*error*/: return 1;
  case 5/*critical*/: return 0;
  }
  return 0;
}

bool SherpaMEs::Initialize(const std::vector<std::string> &args)
{
  p_sherpa = new Sherpa(1);
  std::vector<char*> argv;
  addParameter(argv,"Sherpa");
  addParameter(argv,"EVENTS=0");
  addParameter(argv,"INIT_ONLY=6");
  int level(SherpaVerbosity(spdlog::get("achilles")->level()));
  addParameter(argv,"OUTPUT="+ToString(level));
  // Set beams and PDFs.
  addParameter(argv,"BEAM_1=2212");
  addParameter(argv,"BEAM_2=11");
  addParameter(argv,"BEAM_ENERGY_1=100");
  addParameter(argv,"BEAM_ENERGY_2=100");
  // addParameter(argv,"PDF_SET=None");
  // addParameter(argv,"PDF_LIBRARY=None");
  addParameter(argv,"ME_SIGNAL_GENERATOR=Comix");
  // addParameter(argv,"MODEL=DarkNeutrinoPortal_Dirac_UFO");
  addParameter(argv,"LEPTONIC_CURRENT_MODE=1");
  addParameter(argv,"ALPHAQED_DEFAULT_SCALE=0");
  addParameter(argv,"1/ALPHAQED(MZ)=137");
  // addParameter(argv,"ACTIVE[23]=0");
  addParameter(argv,"ACTIVE[25]=0");
  // addParameter(argv,"ACTIVE[9000005]=0");
  // addParameter(argv,"UFO_PARAM_CARD=parameters.dat");
  addParameter(argv,"PRIORITY[2112]=99");
  addParameter(argv,"PRIORITY[2212]=99");
  addParameter(argv,"MASSIVE[13]=1");
  addParameter(argv,"MASSIVE[15]=1");
  // TODO: Make this something that is passed in
  addParameter(argv,"PRIORITY[1000060120]=99");
  // add additional commandline parameters
  for (const auto &arg: args) addParameter(argv,arg);
  // Initialise Sherpa and return.
  p_sherpa->InitializeTheRun(argv.size(),&argv[0]);
  m_pmap[nlo_type::lo] = new StringProcess_Map();
  RegisterParticles();

  // Clean up memory
  for(auto & arg : argv) delete arg;
  argv.resize(0);
  return true;
}

Process_Base* SherpaMEs::getProcess(Cluster_Amplitude* const ampl) {
  std::string megen = "Comix";
  StringProcess_Map* pm(m_pmap[nlo_type::lo]);
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/"+megen+"/");
  PHASIC::Process_Info pi;
  // pi.m_maxcpl[0]   = pi.m_mincpl[0] = ampl->OrderQCD();
  // pi.m_maxcpl[1]   = pi.m_mincpl[1] = ampl->OrderEW();
  pi.m_megenerator = megen;
  // pi.m_ntchan = 3;
  // pi.m_mtchan = 3;
  pi.m_gpath = "Diagrams";
  for (size_t i(0); i<ampl->NIn(); ++i) {
    Flavour fl(ampl->Leg(i)->Flav().Bar());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_ii.m_ps.emplace_back(fl,"","");
  }
  for (size_t i(ampl->NIn()); i<ampl->Legs().size(); ++i) {
    Flavour fl(ampl->Leg(i)->Flav());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_fi.m_ps.emplace_back(fl,"","");
  }
  msg_Info()<<"SherpaMEs::getProcess: Initializing process ";
  Process_Base* proc = p_sherpa->GetInitHandler()->
    GetMatrixElementHandler()->Generators()->InitializeProcess(pi,false);
  msg_Info()<<" done."<<std::endl;
  if (proc == nullptr) {
    My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
      +"/Process/"+megen+"/");
    My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
      +"/Process/Sherpa/");
    return nullptr;
  }
  m_procs.push_back(proc);
  proc->SetSelector(Selector_Key(nullptr,new Data_Reader(),true));
  proc->SetScale(Scale_Setter_Arguments(MODEL::s_model,
    "VAR{100}{100}", "Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  msg_Info()<<"SherpaMEs::getProcess: Performing tests ";
  if (proc->Get<COMIX::Process_Base>())
    proc->Get<COMIX::Process_Base>()->Tests();
  msg_Info()<<" done."<<std::endl;
  proc->FillProcessMap(&m_pmap);
  proc->ConstructColorMatrix();
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/"+megen+"/");
  My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")
    +"/Process/Sherpa/");
  return proc;
}

bool SherpaMEs::InitializeProcess(const Process_Info &info)
{
  Cluster_Amplitude* ampl = Cluster_Amplitude::New();
  int nqcd(0), nIn(0), cmin(std::numeric_limits<int>::max()), cmax(0);
  for (const auto &initial : info.m_states.begin()->first) {
    ampl->CreateLeg(Vec4D(),Flavour(initial).Bar());
  }
  for (size_t i(0);i<info.m_ids.size();++i) {
    ampl->CreateLeg(Vec4D(),i<1?Flavour(info.m_ids[i]).Bar():
		    Flavour(info.m_ids[i]));
  } 
  for (const auto &final : info.m_states.begin()->second) {
    ampl->CreateLeg(Vec4D(),Flavour(final));
  }
  ampl->SetNIn(2);
  ampl->SetOrderQCD(0);
  ampl->SetOrderEW(info.m_ids.size()-2);
  // TODO: Need to fix this later
  Process_Base::SortFlavours(ampl);
  StringProcess_Map *pm(m_pmap[nlo_type::lo]);
  std::string name(Process_Base::GenerateName(ampl));
  if (pm->find(name)==pm->end()) getProcess(ampl);
  Process_Base* proc(pm->find(name)->second);
  if (proc==nullptr) return false;
  return true;
}

std::map<size_t, long> SherpaMEs::MomentumMap(const std::vector<long> &_fl) const {
    Cluster_Amplitude *ampl(Cluster_Amplitude::New());
    for (size_t i(0);i<_fl.size();++i) {
      Flavour fl(Flavour((long int)(_fl[i])));
      ampl->CreateLeg({},i<2?fl.Bar():fl,ColorID(0,0));
    }
    ampl->SetNIn(2);
    ampl->SetOrderQCD(0);
    ampl->SetOrderEW(ampl->Legs().size()-2);
    Process_Base::SortFlavours(ampl);
    StringProcess_Map *pm(m_pmap[nlo_type::lo]);
    std::string name(Process_Base::GenerateName(ampl));
    spdlog::info("Looking for process");
    if (pm->find(name)==pm->end())
      THROW(fatal_error,"Process not found: "+name);
    Process_Base *proc(pm->find(name)->second);
    auto *singleProcess = proc->Get<COMIX::Single_Process>();
    std::map<size_t, long> mom_map;
    size_t idx = 0;
    for(const auto &flav : singleProcess -> Flavours())
        mom_map[idx++] = flav;

    return mom_map;
}

std::vector<std::unique_ptr<PHASIC::Channels>> SherpaMEs::GenerateChannels(const std::vector<long> &_fl) const {
    Cluster_Amplitude *ampl(Cluster_Amplitude::New());
    for (size_t i(0);i<_fl.size();++i) {
      Flavour fl(Flavour((long int)(_fl[i])));
      ampl->CreateLeg({},i<2?fl.Bar():fl,ColorID(0,0));
    }
    ampl->SetNIn(2);
    ampl->SetOrderQCD(0);
    ampl->SetOrderEW(ampl->Legs().size()-2);
    Process_Base::SortFlavours(ampl);
    StringProcess_Map *pm(m_pmap[nlo_type::lo]);
    std::string name(Process_Base::GenerateName(ampl));
    spdlog::info("Looking for process");
    if (pm->find(name)==pm->end())
      THROW(fatal_error,"Process not found: "+name);
    Process_Base *proc(pm->find(name)->second);
    auto *singleProcess = proc->Get<COMIX::Single_Process>();
    spdlog::info("Process found");

    // Setup external particle channel components
    std::map<int, std::vector<std::shared_ptr<ChannelNode>>> channelComponents;
    std::vector<double> s;
    const auto flavs = singleProcess -> Flavours();
    for(size_t i = 0; i < flavs.size(); ++i) {
        auto node = std::make_shared<ChannelNode>();
        node -> m_pid = flavs[i];
        node -> m_idx = 1 << i;
        channelComponents[(1 << i)].push_back(node);
        s.push_back(sqr(Flavour((kf_code)(flavs[i])).Mass(false)));
    }

    for(unsigned int nset = 2; nset < flavs.size(); ++nset) {
      unsigned int cur = (1 << nset) - 1;
      // Combine all currents
      while(cur < (1 << (flavs.size()))) {
        auto set = SetBits(cur, flavs.size());
        for(unsigned int iset = 1; iset < nset; ++iset) {
          unsigned int idx = (1 << iset) - 1;
          while(idx < (1 << (nset - 1))) {
            unsigned int subCur1 = 0;
            for(unsigned int i = 0; i < flavs.size(); ++i)
              subCur1 += set[i]*((idx >> i) & 1);
            auto subCur2 = cur ^ subCur1;
            // Skip over initial nucleon
            if(SetBit(subCur1, 0) || SetBit(subCur2, 0)) break;
            if(singleProcess -> Combinable(subCur1, subCur2)) {
              // Create new channel component
              for(const auto &subChan1 : channelComponents[subCur1]) {
                for(const auto &subChan2 : channelComponents[subCur2]) {
                  if(cur == (1 << flavs.size()) - 2) {
                    auto node = std::make_shared<ChannelNode>();
                    node -> m_left = subChan1;
                    node -> m_right = subChan2;
                    node -> m_pid = flavs[0];
                    node -> m_idx = cur;
                    channelComponents[cur].push_back(node);
                  } else {
                    for(const auto &elm : singleProcess -> CombinedFlavour(cur)) {
                      auto node = std::make_shared<ChannelNode>();
                      node -> m_left = subChan1;
                      node -> m_right = subChan2;
                      node -> m_pid = static_cast<int>(elm.Kfcode());
                      node -> m_idx = cur;
                      channelComponents[cur].push_back(node);
                    }
                  }
                }
              }

            }
            idx = NextPermutation(idx);
          }
        }
        cur = NextPermutation(cur);
      }
    }

    size_t lid = (1 << _fl.size()) - 2, rid = 2;
    if(channelComponents.find(lid) == channelComponents.end()) throw;
    std::vector<std::unique_ptr<PHASIC::Channels>> channels;
    for(const auto &cur : channelComponents[lid]) {
      auto channel = std::make_unique<GenChannel>(_fl.size(), s);
      channel -> InitializeChannel(cur);
      spdlog::debug("{}", channel -> ToString());
      channels.push_back(std::move(channel));
    }

    return channels;
}

achilles::SherpaMEs::LeptonCurrents SherpaMEs::Calc
(const std::vector<int> &_fl,
 const std::vector<std::array<double, 4> > &p,
 const double &mu2) const
{
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (size_t i(0);i<_fl.size();++i) {
    Vec4D cp(p[i][0],p[i][1],p[i][2],p[i][3]);
    Flavour fl(Flavour((long int)(_fl[i])));
    ampl->CreateLeg(i<2?-cp:cp,i<2?fl.Bar():fl,ColorID(0,0));
  }
  ampl->SetNIn(2);
  ampl->SetOrderQCD(0);
  ampl->SetOrderEW(ampl->Legs().size()-2);
  ampl->SetMuQ2(mu2);
  ampl->SetMuF2(mu2);
  ampl->SetMuR2(mu2);
  Process_Base::SortFlavours(ampl);
  StringProcess_Map *pm(m_pmap[nlo_type::lo]);
  std::string name(Process_Base::GenerateName(ampl));
  if (pm->find(name)==pm->end())
    THROW(fatal_error,"Process not found: "+name);
  Process_Base *proc(pm->find(name)->second);
  // return differntial xs for now
  double res(proc->Differential(*ampl,1|2|4));
  ampl->Delete();

  auto *singleProcess = proc->Get<COMIX::Single_Process>();
  return singleProcess -> LeptonicCurrent();
}

std::vector<FormFactorInfo> achilles::SherpaMEs::FormFactors(int hpid, int vpid) const {
    const std::vector<MODEL::Single_Vertex> &vertices(MODEL::s_model->OriginalVertices());
    std::vector<FormFactorInfo> form_factors;
    for(const auto &vertex : vertices) {
        if(std::find(vertex.in.begin(), vertex.in.end(), -hpid) != vertex.in.end()
                && std::find(vertex.in.begin(), vertex.in.end(), vpid) != vertex.in.end()) {
            for(size_t i = 0; i < vertex.FormFactor.size(); ++i) {
                std::string ff = vertex.FormFactor[i];
                spdlog::trace("For vertex {},{}: found form factor: {}", hpid, vpid, ff);
                std::complex<double> coupling = vertex.Coupling(i);
                FormFactorInfo::Type type;
                if(ff == "F1p") type = FormFactorInfo::Type::F1p;
                else if(ff == "F1n") type = FormFactorInfo::Type::F1n;
                else if(ff == "F2p") type = FormFactorInfo::Type::F2p;
                else if(ff == "F2n") type = FormFactorInfo::Type::F2n;
                else if(ff == "FA") type = FormFactorInfo::Type::FA;
                else if(ff == "FCoh") type = FormFactorInfo::Type::FCoh;
                form_factors.push_back({type, coupling});
            }
        }
    }
    return form_factors;
}

void achilles::SherpaMEs::RegisterParticles() const {
    for(const auto &particleEntry : ATOOLS::s_kftable) {
        auto pid = particleEntry.first;
        auto particle = particleEntry.second;
        static constexpr double to_MeV = 1000;
        const auto mass = particle -> m_mass*to_MeV;
        const auto width = particle -> m_width*to_MeV;
        auto entry = std::make_shared<ParticleInfoEntry>(pid, mass, width, particle->m_icharge,
                                                         particle->m_strong, particle->m_spin,
                                                         particle->m_stable, particle->m_majorana,
                                                         particle->m_massive, particle->m_hadron,
                                                         particle->m_idname, particle->m_antiname);
        achilles::ParticleInfo::Database()[pid] = entry;
    }

    achilles::Database::PrintParticle();
}
