// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0 and gamma spectra at 29 GeV
  class TPC_1985_I205868 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TPC_1985_I205868);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histPhoton, 1, 1, 1);
      book(_histPi    , 2, 1, 1);
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
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");

      for( const Particle& p : ufs.particles()) {
        const int id = p.abspid();
        double xE = p.E()/meanBeamMom;
        switch (id) {
        case 22: // Photons
          _histPhoton->fill(xE);
          break;
        case 111: // Neutral pions
          _histPi->fill(xE);
          break;
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histPhoton, 1./sumOfWeights());
      scale(_histPi    , 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _histPhoton;
    Histo1DPtr _histPi    ;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TPC_1985_I205868);


}
