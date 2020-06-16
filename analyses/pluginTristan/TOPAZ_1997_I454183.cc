// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief N charged vs thrust
  class TOPAZ_1997_I454183 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TOPAZ_1997_I454183);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      ChargedFinalState cfs;
      declare(cfs        , "CFS");
      declare(Thrust(cfs), "Thrust");

      // Book histograms
      book(_p_charged ,3,1,1);
      book(_c_ncharged,1,1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 5 charged FS particles
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 5) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _c_ncharged->fill(cfs.particles().size());
      _p_charged->fill(-log(1.-thrust.thrust()),cfs.particles().size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c_ncharged,1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Profile1DPtr _p_charged;
    CounterPtr _c_ncharged;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TOPAZ_1997_I454183);


}
