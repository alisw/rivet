// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief OPAL strange baryon paper
  /// @author Peter Richardson
  class OPAL_1997_I421977 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_1997_I421977);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_plus,  1, 1, 1);
      book(_h_minus, 2, 1, 1);
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

      for (const Particle& p : ufs.particles()) {
        const int id = p.abspid();
        double xE = p.E()/meanBeamMom;
        if(id==3222)      _h_plus ->fill(xE);
        else if(id==3112) _h_minus->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double fact=1./sumOfWeights();
      scale(_h_plus , fact);
      scale(_h_minus, fact);
      
    }

    //@}


    /// @name Histograms
    //@{
      Histo1DPtr _h_plus,_h_minus;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1997_I421977);


}
