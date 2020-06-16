// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class HRS_1992_I339573 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(HRS_1992_I339573);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_lambda, 1 , 1, 1);
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

      // Final state to get particle spectra
      for( const Particle& p : apply<UnstableFinalState>(event, "UFS").particles(Cuts::abspid==3122)) {
	double xE = p.E()/meanBeamMom;
	_h_lambda->fill(xE);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact = sqr(sqrtS())/GeV2*crossSection()/nanobarn/sumOfWeights();
      scale(_h_lambda, fact); 
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_lambda;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(HRS_1992_I339573);


}
