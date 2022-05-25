// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TriggerUA5.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// UA5 charged multiplicity measurements at 546 GeV
  class UA5_1987_S1640666 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(UA5_1987_S1640666);


    /// Book histograms and initialise projections before the run
    void init() {
      declare(TriggerUA5(), "Trigger");
      declare(ChargedFinalState((Cuts::etaIn(-5.0, 5.0))), "CFS");

      book(_hist_mean_nch ,1, 1, 1);
      book(_hist_nch      ,3, 1, 1);
      book(_sumWPassed, "SumW");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = apply<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;

      _sumWPassed->fill();

      // Count final state particles in several eta regions
      const int Nch = apply<ChargedFinalState>(event, "CFS").size();

      // Fill histograms
      _hist_nch->fill(Nch);
      _hist_mean_nch->fill(_hist_mean_nch->bin(0).xMid(), Nch);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_hist_nch, 1.0 / *_sumWPassed);
      scale(_hist_mean_nch, 1.0 / *_sumWPassed);

    }


  private:

    CounterPtr _sumWPassed;

    Histo1DPtr _hist_mean_nch;
    Histo1DPtr _hist_nch;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(UA5_1987_S1640666, UA5_1987_I244829);

}
