// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/contrib/SoftDrop.hh"

namespace Rivet {



  /// @brief MC validation analysis for jet events
  class MC_JETS : public MC_JetAnalysis {
  public:

    MC_JETS()
      : MC_JetAnalysis("MC_JETS", 4, "Jets")
    {    }


    void init() {
      FinalState fs;
      FastJets jetpro(fs, FastJets::ANTIKT, 0.4);
      const string groomopt = getOption("GROOM", "");
      if (groomopt == "SD") {
        jetpro.addTrf(new fastjet::contrib::SoftDrop(0.0, 0.1));
      } else if (groomopt == "TRIM") {
        jetpro.addTrf(new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
      } else if (groomopt != "") {
        MSG_WARNING("Unknown GROOM=" + groomopt + " option. Not applying jet grooming");
      }
      declare(jetpro, "Jets");
      MC_JetAnalysis::init();
    }


    void analyze(const Event& event) {
      MC_JetAnalysis::analyze(event);
    }


    void finalize() {
      MC_JetAnalysis::finalize();
    }

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_JETS);

}
