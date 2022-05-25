// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// STAR di-hadron correlations in d-Au at 200 GeV
  class STAR_2008_S7993412 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2008_S7993412);


    /// @name Analysis methods
    /// @{

    /// Book projections and histograms
    void init() {
      ChargedFinalState fs((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  1.0*GeV));
      declare(fs, "FS");

      book(_h_Y_jet_trigger ,1, 1, 1);
      book(_h_Y_jet_associated ,2, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Skip if the event is empty
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.empty()) {
        MSG_DEBUG("Skipping event " << numEvents() << " because no final state found ");
        vetoEvent;
      }

      for (const Particle& tp : fs.particles()) {
        const double triggerpT = tp.pT();
        if (triggerpT >= 2.0 && triggerpT < 5.0) {
          int n_associated = 0;
          for (const Particle& ap : fs.particles()) {
            if (!inRange(ap.pT()/GeV, 1.5, triggerpT)) continue;
            if (deltaPhi(tp.phi(), ap.phi()) > 1) continue;
            if (fabs(tp.eta() - ap.eta()) > 1.75) continue;
            n_associated += 1;
          }
          //const double dPhidEta = 2 * 2*1.75;
          //_h_Y_jet_trigger->fill(triggerpT, n_associated/dPhidEta);
          _h_Y_jet_trigger->fill(triggerpT, n_associated);
        }
      }
    }


    /// Finalize
    // void finalize() {    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Profile1DPtr _h_Y_jet_trigger;
    Profile1DPtr _h_Y_jet_associated;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(STAR_2008_S7993412, STAR_2008_I810030);

}
