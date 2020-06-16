// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// NA22 min bias multiplicity distributions.
  class NA22_1986_I18431 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(NA22_1986_I18431);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(ChargedFinalState(), "CFS");

      // Figure out beam type
      const ParticlePair& beam = beams();
      int btype = 0;
      if (beam.first.pid() == PID::PIPLUS && beam.second.pid() == PID::PROTON)
        btype = 1;
      else if (beam.first.pid() == PID::KPLUS && beam.second.pid() == PID::PROTON)
        btype = 2;
      else if (beam.first.pid() == PID::PROTON && beam.second.pid() == PID::PROTON)
        btype = 3;
      else {
        MSG_ERROR("Beam error: Not compatible!");
        return;
      }

      // Book histo for appropriate beam type
      book(_h_mult, btype, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      size_t nfs = apply<ChargedFinalState>(event, "CFS").size();
      if (nfs >= 30)  nfs = 30;
      _h_mult->fill(nfs);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mult, 2.0); // normalize to 2 to account for bin width
    }

    ///@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_mult;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(NA22_1986_I18431);


}
