// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief Demonstrate use of the Cutflow(s) classes
  class EXAMPLE_CUTFLOW : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EXAMPLE_CUTFLOW);


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      // Projections
      const DirectFinalState leps(Cuts::abseta < 4 && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      declare(leps, "Leptons");

      const FinalState fs(Cuts::abseta < 4);
      const FastJets jets(fs, FastJets::Algo::ANTIKT, 0.4);
      declare(jets, "Jets");

      // Histograms
      // book(_s_cutflow, "cutflow");

      // Cut-flows
      _cutflows.addCutflow("JetLep", {"Jets", "Nlep", "pTlep1", "yLep1"});
      _cutflows.addCutflow("DiBjet", {"Jets", "0Lep", "2Jet", "Nbjets", });

    }


    /// Do the analysis
    void analyze(const Event& event) {

      _cutflows.fillinit();

      const Jets jets = apply<JetFinder>(event, "Jets").jetsByPt();
      const Particles leps = apply<FinalState>(event, "Leptons").particlesByPt();

      if (jets.empty()) vetoEvent;
      _cutflows.fill(1);

      if (!leps.empty()) {
        _cutflows["JetLep"].fill(2);
        if (!_cutflows["JetLep"].fillnext(leps[0].pT() < 20*GeV)) vetoEvent;
        _cutflows["JetLep"].fillnext(leps[0].absrap() < 2.5);
      } else {
        Cutflow& cf = _cutflows["DiBjet"];
        cf.fillnext();
        const Jets bjets = filter_select(jets, hasBTag(Cuts::pT > 5*GeV));
        cf.fillnext({jets.size() >= 2, bjets.size() == 2});
      }

    }

    /// Finalize
    void finalize() {
      MSG_INFO("Cut-flow:\n" << _cutflows);
    }

    //@}


    /// Cut-flow counters
    Cutflows _cutflows;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EXAMPLE_CUTFLOW);

}
