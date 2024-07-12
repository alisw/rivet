// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief OPAL charged particle fragmentation functions
  ///
  /// @author Peter Richardson
  class OPAL_1994_S2927284 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1994_S2927284);


    /// @name Analysis methods
    /// @{

    void analyze(const Event& e) {

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& fs = apply<FinalState>(e, "FS");
      if (fs.particles().size() < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      for (const Particle& p : fs.particles()) {
        int id = p.abspid();
        // charged pions
        if (id == PID::PIPLUS) {
          _histXpPiPlus->fill(p.p3().mod());
        } else if(id == PID::KPLUS) {
          _histXpKPlus->fill(p.p3().mod());
        } else if(id == PID::PROTON) {
          _histXpProton->fill(p.p3().mod());
        }
      }
    }


    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      book(_histXpPiPlus , 1, 1, 1);
      book(_histXpKPlus  , 2, 1, 1);
      book(_histXpProton , 3, 1, 1);
    }


    /// Finalize
    void finalize() {
      scale(_histXpPiPlus,1./sumOfWeights());
      scale(_histXpKPlus ,1./sumOfWeights());
      scale(_histXpProton,1./sumOfWeights());
    }

    /// @}


  private:

    Histo1DPtr _histXpPiPlus;
    Histo1DPtr _histXpKPlus;
    Histo1DPtr _histXpProton;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(OPAL_1994_S2927284, OPAL_1994_I372772);

}
