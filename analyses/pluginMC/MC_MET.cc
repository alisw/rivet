// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {



  /// @brief MC validation analysis for truth-MET measurement
  /// @todo Add plots for MET based on prompt invisibles
  class MC_MET : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(MC_MET);


    void init() {
      const FinalState inclfs;
      FinalState calofs(Cuts::abseta < 5);
      declare(MissingMomentum(inclfs), "InclMET");
      declare(MissingMomentum(calofs), "CaloMET");

      declare(InvisibleFinalState(), "InvisibleFS");
      declare(InvisibleFinalState(true), "PromptInvisibleFS");

      book(_h["met_incl"], "met_incl", logspace(50, 10, sqrtS()/GeV/5));
      book(_h["met_calo"], "met_calo", logspace(50, 10, sqrtS()/GeV/5));
      book(_h["set_incl"], "set_incl", logspace(50, 10, sqrtS()/GeV/3));
      book(_h["set_calo"], "set_calo", logspace(50, 10, sqrtS()/GeV/3));

      book(_h["pT_inv"],   "pT_inv",   logspace(50, 10, sqrtS()/GeV/5));
      book(_h["mass_inv"], "mass_inv", logspace(100, 10, sqrtS()/GeV/5));
      book(_h["rap_inv"],  "rap_inv",   50, -5., 5.);

      book(_h["pT_promptinv"],   "pT_promptinv",   logspace(50, 10, sqrtS()/GeV/5));
      book(_h["mass_promptinv"], "mass_promptinv", logspace(100, 10, sqrtS()/GeV/5));
      book(_h["rap_promptinv"],  "rap_promptinv",  50, -5., 5.);
    }


    void analyze(const Event& event) {

      const MissingMomentum& mmincl = apply<MissingMomentum>(event, "InclMET");
      _h["met_incl"]->fill(mmincl.met()/GeV);
      _h["set_incl"]->fill(mmincl.set()/GeV);

      const MissingMomentum& mmcalo = apply<MissingMomentum>(event, "CaloMET");
      _h["met_calo"]->fill(mmcalo.met()/GeV);
      _h["set_calo"]->fill(mmcalo.set()/GeV);

      // Get the invisible final state particles
      const Particles& invisibles = apply<InvisibleFinalState>(event, "InvisibleFS").particlesByPt();
      const Particles& promptinvisibles = apply<InvisibleFinalState>(event, "PromptInvisibleFS").particlesByPt();

      FourMomentum invsum;
      for (const Particle& p : invisibles) {
        invsum += p.momentum();
      }
      _h["pT_inv"]->fill(invsum.pT()/GeV);
      _h["mass_inv"]->fill(invsum.mass()/GeV);
      _h["rap_inv"]->fill(invsum.rapidity());

      FourMomentum promptinvsum;
      for (const Particle& p : promptinvisibles) {
        promptinvsum += p.momentum();
      }
      _h["pT_promptinv"]->fill(promptinvsum.pT()/GeV);
      _h["mass_promptinv"]->fill(promptinvsum.mass()/GeV);
      _h["rap_promptinv"]->fill(promptinvsum.rapidity());

    }


    void finalize() {
      const double sf = crossSectionPerEvent()/picobarn;
      scale(_h, sf);
    }


  private:

    map<string, Histo1DPtr> _h;

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_MET);
}
