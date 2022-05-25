// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/TriggerCDFRun0Run1.hh"

namespace Rivet {


  /// @brief CDF track \f$ p_\perp \f$ distributions at 630 and 1800 GeV
  class CDF_1988_S1865951 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_1988_S1865951);


    /// @name Analysis methods
    //@{

    /// Book histograms and set up projections
    void init() {
      // Set up projections
      declare(TriggerCDFRun0Run1(), "Trigger");
      const ChargedFinalState cfs((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  0.4*GeV));
      declare(cfs, "CFS");

      // Book histo
      if (isCompatibleWithSqrtS(1800)) {
        book(_hist_pt ,1, 1, 1);
      } else if (isCompatibleWithSqrtS(630)) {
        book(_hist_pt ,2, 1, 1);
      }

      book(_sumWTrig, "sumWTrig");

    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Trigger
      const bool trigger = apply<TriggerCDFRun0Run1>(event, "Trigger").minBiasDecision();
      if (!trigger) vetoEvent;
      _sumWTrig->fill();

      const FinalState& trackfs = apply<ChargedFinalState>(event, "CFS");
      for (Particle p : trackfs.particles()) {
        const double pt = p.pT()/GeV;
        // Effective weight for d3sig/dp3 = weight / ( Delta eta * 2pi * pt ), with Delta(eta) = 2
        const double eff_weight = 1.0/(2*2*TWOPI*pt);
        _hist_pt->fill(pt, eff_weight);
      }
    }


    /// Scale histos
    void finalize() {
      scale(_hist_pt, crossSectionPerEvent()/millibarn);
    }

    //@}


  private:

    /// Counter
    CounterPtr _sumWTrig;

    /// Histo
    Histo1DPtr _hist_pt;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_1988_S1865951, CDF_1988_I263320);

}
