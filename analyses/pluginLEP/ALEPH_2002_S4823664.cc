// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ALEPH eta/omega fragmentation function paper
  ///
  /// @author Peter Richardson
  class ALEPH_2002_S4823664 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_2002_S4823664);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXpEta   , 2, 1, 2);
      book(_histXpOmega , 3, 1, 2);
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
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");

      for (const Particle& p : ufs.particles()) {
        if(p.abspid()==221) {
          double xp = p.p3().mod()/meanBeamMom;
          _histXpEta->fill(xp);
        }
        else if(p.abspid()==223) {
          double xp = p.p3().mod()/meanBeamMom;
          _histXpOmega->fill(xp);
        }
      }
    }


    void finalize() {
      scale(_histXpEta  , 1./sumOfWeights());
      scale(_histXpOmega, 1./sumOfWeights());
    }

    /// @}


  private:

    Histo1DPtr _histXpEta;
    Histo1DPtr _histXpOmega;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(ALEPH_2002_S4823664, ALEPH_2002_I569165);

}
