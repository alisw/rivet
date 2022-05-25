// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief DELPHI rho,f_0 and f_2 fragmentation function paper
  ///
  /// @author Peter Richardson
  class DELPHI_1999_S3960137 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1999_S3960137);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXpRho , 1, 1, 1);
      book(_histXpf0  , 1, 1, 2);
      book(_histXpf2  , 1, 1, 3);
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
        const int id = p.abspid();
        double xp = p.p3().mod()/meanBeamMom;
        switch (id) {
        case 9010221:
          _histXpf0->fill(xp);
          break;
        case 225:
          _histXpf2->fill(xp);
          break;
        case 113:
          _histXpRho->fill(xp);
          break;
        }
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXpf0 , 1./sumOfWeights());
      scale(_histXpf2 , 1./sumOfWeights());
      scale(_histXpRho, 1./sumOfWeights());
    }

    /// @}


  private:

      Histo1DPtr _histXpf0;
      Histo1DPtr _histXpf2;
      Histo1DPtr _histXpRho;
  };



  RIVET_DECLARE_ALIASED_PLUGIN(DELPHI_1999_S3960137, DELPHI_1999_I482816);

}
