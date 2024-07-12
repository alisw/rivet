// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief OPAL strange baryon paper
  ///
  /// @author Peter Richardson
  class OPAL_1997_S3396100 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1997_S3396100);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_histXpLambda         , 1, 1, 1);
      book(_histXiLambda         , 2, 1, 1);
      book(_histXpXiMinus        , 3, 1, 1);
      book(_histXiXiMinus        , 4, 1, 1);
      book(_histXpSigma1385Plus  , 5, 1, 1);
      book(_histXiSigma1385Plus  , 6, 1, 1);
      book(_histXpSigma1385Minus , 7, 1, 1);
      book(_histXiSigma1385Minus , 8, 1, 1);
      book(_histXpXi1530         , 9, 1, 1);
      book(_histXiXi1530         ,10, 1, 1);
      book(_histXpLambda1520     ,11, 1, 1);
      book(_histXiLambda1520     ,12, 1, 1);
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
        if (!inRange(id, 3000, 3999)) continue;

        const double xE = p.E()/meanBeamMom;
        const double xp = p.p3().mod()/(2*meanBeamMom);
        const double xi = -log(xp);

        switch (id) {
        case 3312:
          _histXpXiMinus->fill(xE);
          _histXiXiMinus->fill(xi);
          break;
        case 3224:
          _histXpSigma1385Plus->fill(xE);
          _histXiSigma1385Plus->fill(xi);
          break;
        case 3114:
          _histXpSigma1385Minus->fill(xE);
          _histXiSigma1385Minus->fill(xi);
          break;
        case 3122:
          _histXpLambda->fill(xE);
          _histXiLambda->fill(xi);
          break;
        case 3324:
          _histXpXi1530->fill(xE);
          _histXiXi1530->fill(xi);
          break;
        case 3124:
          _histXpLambda1520->fill(xE);
          _histXiLambda1520->fill(xi);
          break;
        }
      }
    }


    /// Finalize
    void finalize() {
      double fact=1./sumOfWeights();
      scale(_histXpLambda        , fact);
      scale(_histXiLambda        , fact);
      scale(_histXpXiMinus       , fact);
      scale(_histXiXiMinus       , fact);
      scale(_histXpSigma1385Plus , fact);
      scale(_histXiSigma1385Plus , fact);
      scale(_histXpSigma1385Minus, fact);
      scale(_histXiSigma1385Minus, fact);
      scale(_histXpXi1530        , fact);
      scale(_histXiXi1530        , fact);
      scale(_histXpLambda1520    , fact);
      scale(_histXiLambda1520    , fact);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _histXpLambda;
    Histo1DPtr _histXiLambda;
    Histo1DPtr _histXpXiMinus;
    Histo1DPtr _histXiXiMinus;
    Histo1DPtr _histXpSigma1385Plus;
    Histo1DPtr _histXiSigma1385Plus;
    Histo1DPtr _histXpSigma1385Minus;
    Histo1DPtr _histXiSigma1385Minus;
    Histo1DPtr _histXpXi1530;
    Histo1DPtr _histXiXi1530;
    Histo1DPtr _histXpLambda1520;
    Histo1DPtr _histXiLambda1520;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(OPAL_1997_S3396100, OPAL_1997_I421978);

}
