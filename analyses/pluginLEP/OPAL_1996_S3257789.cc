// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief OPAL J/Psi fragmentation function paper
  /// @author Peter Richardson
  class OPAL_1996_S3257789 : public Analysis {
  public:

    /// Constructor
    OPAL_1996_S3257789()
      : Analysis("OPAL_1996_S3257789")
    {}


    /// @name Analysis methods
    //@{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXpJPsi   , 1, 1, 1);
      book(_multJPsi     , 2, 1, 1);
      book(_multPsiPrime , 2, 1, 2);
      book(_weightSum,"_weightSum");
    }


    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");

      for (const Particle& p : ufs.particles()) {
        if(p.abspid()==443) {
          double xp = p.p3().mod()/meanBeamMom;
          _histXpJPsi->fill(xp);
          _multJPsi->fill(91.2);
          _weightSum->fill();
        }
        else if(p.abspid()==100443) {
          _multPsiPrime->fill(91.2);
        }
      }
    }


    /// Finalize
    void finalize() {
      if(_weightSum->val()>0.)
        scale(_histXpJPsi  , 0.1/ *_weightSum);
      scale(_multJPsi    , 1./sumOfWeights());
      scale(_multPsiPrime, 1./sumOfWeights());
    }

    //@}


  private:

    CounterPtr _weightSum;
    Histo1DPtr _histXpJPsi;
    Histo1DPtr _multJPsi;
    Histo1DPtr _multPsiPrime;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1996_S3257789);

}
