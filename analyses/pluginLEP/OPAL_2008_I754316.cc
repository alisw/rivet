// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/GammaGammaFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class OPAL_2008_I754316 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2008_I754316);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // get the hadronic final state
      const GammaGammaKinematics& gammakin = declare(GammaGammaKinematics(), "Kinematics");
      const FinalState & fs = declare(GammaGammaFinalState(gammakin), "FS");
      declare(FastJets(fs, FastJets::KT,1.),"Jets");

      // Book histograms
      book(_h_y1,1, 1, 1);
      book(_h_y2,2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 5*GeV and Cuts::abseta < 1.5);
      if(jets.empty()) vetoEvent;
      for(const Jet & jet : jets) {
      	_h_y2->fill(jet.pT());
      	if(abs(jet.eta())<1.0)
      	  _h_y1->fill(jet.pT());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_y1, crossSection()/picobarn/sumOfWeights());
      scale(_h_y2, crossSection()/picobarn/sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_y1, _h_y2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2008_I754316);


}
