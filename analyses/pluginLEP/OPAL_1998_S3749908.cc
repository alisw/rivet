// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief OPAL photon/light meson paper
  ///
  /// @author Peter Richardson
  class OPAL_1998_S3749908 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1998_S3749908);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXePhoton   , 2, 1, 1);
      book(_histXiPhoton   , 3, 1, 1);
      book(_histXePi       , 4, 1, 1);
      book(_histXiPi       , 5, 1, 1);
      book(_histXeEta      , 6, 1, 1);
      book(_histXiEta      , 7, 1, 1);
      book(_histXeRho      , 8, 1, 1);
      book(_histXiRho      , 9, 1, 1);
      book(_histXeOmega    ,10, 1, 1);
      book(_histXiOmega    ,11, 1, 1);
      book(_histXeEtaPrime ,12, 1, 1);
      book(_histXiEtaPrime ,13, 1, 1);
      book(_histXeA0       ,14, 1, 1);
      book(_histXiA0       ,15, 1, 1);
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
        double xi = -log(p.p3().mod()/meanBeamMom);
        double xE = p.E()/meanBeamMom;
        switch (id) {
        case 22: // Photons
          _histXePhoton->fill(xE);
          _histXiPhoton->fill(xi);
          break;
        case 111: // Neutral pions
          _histXePi->fill(xE);
          _histXiPi->fill(xi);
          break;
        case 221: // eta
          _histXeEta->fill(xE);
          _histXiEta->fill(xi);
          break;
        case 213: // Charged rho (770)
          _histXeRho->fill(xE);
          _histXiRho->fill(xi);
          break;
        case 223: // omega (782)
          _histXeOmega->fill(xE);
          _histXiOmega->fill(xi);
          break;
        case 331: // eta' (958)
          _histXeEtaPrime->fill(xE);
          _histXiEtaPrime->fill(xi);
          break;
        case 9000211: // Charged a_0 (980)
          _histXeA0->fill(xE);
          _histXiA0->fill(xi);
          break;
        }
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXePhoton  , 1./sumOfWeights());
      scale(_histXiPhoton  , 1./sumOfWeights());
      scale(_histXePi      , 1./sumOfWeights());
      scale(_histXiPi      , 1./sumOfWeights());
      scale(_histXeEta     , 1./sumOfWeights());
      scale(_histXiEta     , 1./sumOfWeights());
      scale(_histXeRho     , 1./sumOfWeights());
      scale(_histXiRho     , 1./sumOfWeights());
      scale(_histXeOmega   , 1./sumOfWeights());
      scale(_histXiOmega   , 1./sumOfWeights());
      scale(_histXeEtaPrime, 1./sumOfWeights());
      scale(_histXiEtaPrime, 1./sumOfWeights());
      scale(_histXeA0      , 1./sumOfWeights());
      scale(_histXiA0      , 1./sumOfWeights());
    }

    /// @}


  private:

    /// @{
    Histo1DPtr _histXePhoton;
    Histo1DPtr _histXiPhoton;
    Histo1DPtr _histXePi;
    Histo1DPtr _histXiPi;
    Histo1DPtr _histXeEta;
    Histo1DPtr _histXiEta;
    Histo1DPtr _histXeRho;
    Histo1DPtr _histXiRho;
    Histo1DPtr _histXeOmega;
    Histo1DPtr _histXiOmega;
    Histo1DPtr _histXeEtaPrime;
    Histo1DPtr _histXiEtaPrime;
    Histo1DPtr _histXeA0;
    Histo1DPtr _histXiA0;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(OPAL_1998_S3749908, OPAL_1998_I470419);

}
