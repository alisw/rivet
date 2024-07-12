// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Jet production in photon-photon collisions at 198 GeV
  class L3_2004_I661114 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_2004_I661114);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // get the hadronic final state
      const GammaGammaKinematics& diskin = declare(GammaGammaKinematics(), "Kinematics");
      const FinalState & fs = declare(GammaGammaFinalState(diskin), "FS");
      declare(FastJets(fs, FastJets::KT,1.),"Jets");

      // Book histograms
      book(_h_y, 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 3*GeV and Cuts::abseta < 1.0);
      if(jets.empty()) vetoEvent;
      for(const Jet & jet : jets) {
      	_h_y->fill(jet.pT());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_y, crossSection()/picobarn/sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_y;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(L3_2004_I661114);


}
