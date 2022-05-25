// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief DELPHI strange baryon paper
  ///
  /// @author Hendrik Hoeth
  class DELPHI_1995_S3137023 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1995_S3137023);


    /// @name Analysis methods
    /// @{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      book(_histXpXiMinus       ,2, 1, 1);
      book(_histXpSigma1385Plus ,3, 1, 1);
      book(_weightedTotalNumXiMinus, "_weightedTotalNumXiMinus");
      book(_weightedTotalNumSigma1385Plus, "_weightedTotalNumSigma1385Plus");

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
        switch (id) {
        case 3312:
          _histXpXiMinus->fill(p.p3().mod()/meanBeamMom);
          _weightedTotalNumXiMinus->fill();
          break;
        case 3114: case 3224:
          _histXpSigma1385Plus->fill(p.p3().mod()/meanBeamMom);
          _weightedTotalNumSigma1385Plus->fill();
          break;
        }
      }

    }


    /// Finalize
    void finalize() {
      normalize(_histXpXiMinus       , dbl(*_weightedTotalNumXiMinus)/sumOfWeights());
      normalize(_histXpSigma1385Plus , dbl(*_weightedTotalNumSigma1385Plus)/sumOfWeights());
    }

    /// @}


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.
    CounterPtr _weightedTotalNumXiMinus;
    CounterPtr _weightedTotalNumSigma1385Plus;

    Histo1DPtr _histXpXiMinus;
    Histo1DPtr _histXpSigma1385Plus;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(DELPHI_1995_S3137023, DELPHI_1995_I394716);

}
