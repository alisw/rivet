// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// @brief Thrust, heavy jet mass, and y3 at 58 GeV
  class TOPAZ_1993_I361661 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TOPAZ_1993_I361661);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs, "FS");
      declare(FastJets(fs, FastJets::DURHAM, 0.7), "DurhamJets");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      book(_h_thrust, 1, 1, 1);
      book(_h_rho   , 2, 1, 1);
      book(_h_y23   , 3, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      _h_thrust->fill(-log(1.-thrust.thrust()));
      // jet mass
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      _h_rho->fill(-log(hemi.scaledM2high()));
      // Jets
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      if(numParticles>=3)
        if (durjet.clusterSeq()) _h_y23->fill(-log(durjet.clusterSeq()->exclusive_ymerge_max(2)));
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_thrust);
      normalize(_h_rho);
      normalize(_h_y23); 

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_thrust,_h_rho,_h_y23;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TOPAZ_1993_I361661);


}
