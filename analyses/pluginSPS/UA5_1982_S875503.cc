// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  /// @brief UA5 multiplicity and \f$ \eta \f$ distributions
  class UA5_1982_S875503 : public Analysis {
  public:

    /// Default constructor
    UA5_1982_S875503() : Analysis("UA5_1982_S875503") {
    }


    /// @name Analysis methods
    //@{

    /// Set up projections and book histos
    void init() {
      declare(TriggerUA5(), "Trigger");
      declare(ChargedFinalState((Cuts::etaIn(-3.5, 3.5))), "CFS");

      // Book histos based on pp or ppbar beams
      if (beamIds().first == beamIds().second) {
        book(_hist_nch ,2,1,1);
        book(_hist_eta ,3,1,1);
      } else {
        book(_hist_nch ,2,1,2);
        book(_hist_eta ,4,1,1);
      }
      book(_sumWTrig, "sumW");

    }


    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = apply<TriggerUA5>(event, "Trigger");
      if (!trigger.nsdDecision()) vetoEvent;
      _sumWTrig->fill();

      // Get tracks
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // Fill mean charged multiplicity histos
      _hist_nch->fill(_hist_nch->bin(0).xMid(), cfs.size());

      // Iterate over all tracks and fill eta histograms
      for (const Particle& p : cfs.particles()) {
        const double eta = p.abseta();
        _hist_eta->fill(eta);
      }

    }


    void finalize() {
      /// @todo Why the factor of 2 on Nch for ppbar?
      if (beamIds().first == beamIds().second) {
        scale(_hist_nch, 1.0 / *_sumWTrig);
      } else {
        scale(_hist_nch, 0.5 / *_sumWTrig);
      }
      scale(_hist_eta, 0.5 / *_sumWTrig);
    }

    //@}


  private:

    /// @name Counters
    //@{
    CounterPtr _sumWTrig;
    //@}

    /// @name Histogram collections
    //@{
    Histo1DPtr _hist_nch;
    Histo1DPtr _hist_eta;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(UA5_1982_S875503);

}
