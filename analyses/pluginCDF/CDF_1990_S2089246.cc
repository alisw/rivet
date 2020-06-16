// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerCDFRun0Run1.hh"

namespace Rivet {


  /// @brief CDF pseudorapidity analysis at 630 and 1800 GeV
  /// @author Andy Buckley
  class CDF_1990_S2089246 : public Analysis {
  public:

    /// Constructor
    CDF_1990_S2089246()
      : Analysis("CDF_1990_S2089246")
    {
    }


    /// @name Analysis methods
    //@{

    void init() {
      // Setup projections
      declare(TriggerCDFRun0Run1(), "Trigger");
      declare(ChargedFinalState((Cuts::etaIn(-3.5, 3.5))), "CFS");

      // Book histo
      if (fuzzyEquals(sqrtS()/GeV, 1800, 1E-3)) {
        book(_hist_eta ,3, 1, 1);
      } else if (fuzzyEquals(sqrtS()/GeV, 630, 1E-3)) {
        book(_hist_eta ,4, 1, 1);
      }
      book(_sumWTrig, "sumWTrig");
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const bool trigger = apply<TriggerCDFRun0Run1>(event, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;
      _sumWTrig->fill();

      // Loop over final state charged particles to fill eta histos
      const FinalState& fs = apply<FinalState>(event, "CFS");
      for (const Particle& p : fs.particles()) {
        const double eta = p.eta();
        _hist_eta->fill(fabs(eta));
      }
    }


    /// Finalize
    void finalize() {
      // Divide through by num events to get d<N>/d(eta) in bins
      // Factor of 1/2 for |eta| -> eta
      scale(_hist_eta, 0.5/ *_sumWTrig);
    }

    //@}


  private:

    /// @name Weight counter
    //@{
    CounterPtr _sumWTrig;
    //@}

    /// @name Histogram collections
    //@{
    Histo1DPtr _hist_eta;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_1990_S2089246);

}
