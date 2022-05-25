// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief K*0, Phi and production
  class DELPHI_1996_I420528 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1996_I420528);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_Kstar,1, 1, 1);
      book(_h_phi  ,3, 1, 1);

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
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles(Cuts::abspid==313 or Cuts::pid==333)) {
        const int id = p.abspid();
        double xp = p.p3().mod()/meanBeamMom;
        switch (id) {
        case 333:
          _h_phi->fill(xp);
          break;
        case 313:
          _h_Kstar->fill(xp);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Kstar , 1./sumOfWeights());
      scale(_h_phi , 1./sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Kstar, _h_phi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_1996_I420528);


}
