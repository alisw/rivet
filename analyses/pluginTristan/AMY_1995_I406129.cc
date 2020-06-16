// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JadePlugin.hh"

namespace fastjet {

class P_scheme : public JetDefinition::Recombiner {
 public:
  std::string description() const {return "";}
  void recombine(const PseudoJet & pa, const PseudoJet & pb,
		 PseudoJet & pab) const {
    PseudoJet tmp = pa + pb;
    double E = sqrt(tmp.px()*tmp.px() + tmp.py()*tmp.py() + tmp.pz()*tmp.pz());
    pab.reset_momentum(tmp.px(), tmp.py(), tmp.pz(), E);
  }
  void preprocess(PseudoJet & p) const {
    double E = sqrt(p.px()*p.px() + p.py()*p.py() + p.pz()*p.pz());
    p.reset_momentum(p.px(), p.py(), p.pz(), E);
  }
  ~P_scheme() { }
};

class E0_scheme : public JetDefinition::Recombiner {
 public:
  std::string description() const {return "";}
  void recombine(const PseudoJet & pa, const PseudoJet & pb,
		 PseudoJet & pab) const {
    PseudoJet tmp = pa + pb;
    double fact = tmp.E()/sqrt(tmp.px()*tmp.px()+tmp.py()*tmp.py()+tmp.pz()*tmp.pz());
    pab.reset_momentum(fact*tmp.px(), fact*tmp.py(), fact*tmp.pz(), tmp.E());
  }
  void preprocess(PseudoJet & p) const {
    double fact = p.E()/sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz());
  
    p.reset_momentum(fact*p.px(), fact*p.py(), fact*p.pz(), p.E());
  }
  ~E0_scheme() { }
};

}


namespace Rivet {


  /// @brief AMY jets at
  class AMY_1995_I406129 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(AMY_1995_I406129);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      FinalState fs;
      declare(fs, "FS");
      // Book histograms
      book(_h_jade_P , 2, 1, 1);
      book(_h_jade_E , 3, 1, 1);
      book(_h_jade_E0, 6, 1, 1);
      book(_h_durham , 4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles particles =  apply<FinalState>(event, "FS").particles();
      MSG_DEBUG("Num particles = " << particles.size());
      PseudoJets pjs;
      double mpi=.13957;
      for (const Particle & p : particles) {
	Vector3 mom = p.p3();
	double energy = p.E();
	if(PID::isCharged(p.pid())) {
	  energy = sqrt(mom.mod2()+sqr(mpi));
	}
	else {
	  double fact = energy/mom.mod();
	  mom *=fact;
	}
	pjs.push_back(fastjet::PseudoJet(mom.x(),mom.y(),mom.z(),energy));
      }
      // durham
      fastjet::JetDefinition durDef(fastjet::ee_kt_algorithm, fastjet::E_scheme);
      fastjet::ClusterSequence durham(pjs,durDef);
      double y_23 = durham.exclusive_ymerge_max(2);
      _h_durham->fill(y_23);
      // jade e-scheme
      fastjet::JetDefinition::Plugin *plugin = new fastjet::JadePlugin();
      fastjet::JetDefinition jadeEDef(plugin);
      jadeEDef.set_recombination_scheme(fastjet::E_scheme);
      fastjet::ClusterSequence jadeE(pjs,jadeEDef);
      y_23 = jadeE.exclusive_ymerge_max(2);
      _h_jade_E->fill(y_23);
      // jade p-scheme
      fastjet::P_scheme p_scheme;
      fastjet::JetDefinition jadePDef(plugin);
      jadePDef.set_recombiner(&p_scheme);
      fastjet::ClusterSequence jadeP(pjs,jadePDef);
      y_23 = jadeP.exclusive_ymerge_max(2);
      _h_jade_P->fill(y_23);
      // jade E0-scheme
      fastjet::E0_scheme e0_scheme;
      fastjet::JetDefinition jadeE0Def(plugin);
      jadeE0Def.set_recombiner(&e0_scheme);
      fastjet::ClusterSequence jadeE0(pjs,jadeE0Def);
      y_23 = jadeE0.exclusive_ymerge_max(2);
      _h_jade_E0->fill(y_23);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_jade_E , 1./sumOfWeights());
      scale(_h_jade_E0, 1./sumOfWeights());
      scale(_h_jade_P , 1./sumOfWeights());
      scale(_h_durham , 1./sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_jade_P, _h_jade_E, _h_jade_E0, _h_durham; 
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(AMY_1995_I406129);


}
