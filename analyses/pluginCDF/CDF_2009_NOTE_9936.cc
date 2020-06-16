// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerCDFRun2.hh"

namespace Rivet {


  class CDF_2009_NOTE_9936 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_2009_NOTE_9936()
      : Analysis("CDF_2009_NOTE_9936")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(TriggerCDFRun2(), "Trigger");

      declare(ChargedFinalState((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  0.4*GeV)), "CFS");

      book(_hist_nch ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // MinBias Trigger
      const bool trigger = apply<TriggerCDFRun2>(event, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;

      // Get events charged multiplicity and fill histogram
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      _hist_nch->fill(cfs.size());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_hist_nch);
    }

    //@}

  private:

    Histo1DPtr _hist_nch;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2009_NOTE_9936);

}
