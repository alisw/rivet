// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  /// UA5 min bias charged multiplicities in central \f$ \eta \f$ ranges
  class UA5_1989_S1926373 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(UA5_1989_S1926373);


    /// @name Analysis methods
    /// @{

    /// Book histograms and projections
    void init() {
      declare(TriggerUA5(), "Trigger");
      declare(ChargedFinalState((Cuts::etaIn(-0.5, 0.5))), "CFS05");
      declare(ChargedFinalState((Cuts::etaIn(-1.5, 1.5))), "CFS15");
      declare(ChargedFinalState((Cuts::etaIn(-3.0, 3.0))), "CFS30");
      declare(ChargedFinalState((Cuts::etaIn(-5.0, 5.0))), "CFS50");

      // NB. _hist_nch and _hist_ncheta50 use the same data but different binning
      if (isCompatibleWithSqrtS(200)) {
        book(_hist_nch        ,1, 1, 1);
        book(_hist_nch_eta05  ,3, 1, 1);
        book(_hist_nch_eta15  ,4, 1, 1);
        book(_hist_nch_eta30  ,5, 1, 1);
        book(_hist_nch_eta50  ,6, 1, 1);
        book(_hist_mean_nch   ,11, 1, 1);
      } else if (isCompatibleWithSqrtS(900)) {
        book(_hist_nch        ,2, 1, 1);
        book(_hist_nch_eta05  ,7, 1, 1);
        book(_hist_nch_eta15  ,8, 1, 1);
        book(_hist_nch_eta30  ,9, 1, 1);
        book(_hist_nch_eta50  ,10, 1, 1);
        book(_hist_mean_nch   ,12, 1, 1);
      }
      book(_sumWPassed, "SumW");
      /// @todo Moments of distributions
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = apply<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;

      _sumWPassed->fill();

      // Count final state particles in several eta regions
      const int numP05 = apply<ChargedFinalState>(event, "CFS05").size();
      const int numP15 = apply<ChargedFinalState>(event, "CFS15").size();
      const int numP30 = apply<ChargedFinalState>(event, "CFS30").size();
      const int numP50 = apply<ChargedFinalState>(event, "CFS50").size();

      // Fill histograms
      _hist_nch->fill(numP50);
      _hist_nch_eta05->fill(numP05);
      _hist_nch_eta15->fill(numP15);
      _hist_nch_eta30->fill(numP30);
      _hist_nch_eta50->fill(numP50);
      _hist_mean_nch->fill(_hist_mean_nch->bin(0).xMid(), numP50);
    }


    void finalize() {
      scale(_hist_nch, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta05, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta15, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta30, 1.0 / *_sumWPassed);
      scale(_hist_nch_eta50, 1.0 / *_sumWPassed);
      scale(_hist_mean_nch, 1.0 / *_sumWPassed);
    }

    /// @}


  private:

    /// Weight counter
    CounterPtr _sumWPassed;

    /// @name Histograms
    /// @{
    Histo1DPtr _hist_nch;
    Histo1DPtr _hist_nch_eta05;
    Histo1DPtr _hist_nch_eta15;
    Histo1DPtr _hist_nch_eta30;
    Histo1DPtr _hist_nch_eta50;
    Histo1DPtr _hist_mean_nch;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(UA5_1989_S1926373, UA5_1989_I267179);

}
