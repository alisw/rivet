// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class OPAL_2003_I595335 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2003_I595335);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const ChargedFinalState cfs = ChargedFinalState();
      declare(cfs, "CFS");
      declare(Thrust(cfs), "Thrust");

      book(_h_pT_in , 1, 1, 1);
      book(_h_pT_out, 2, 1, 1);
      book(_h_y     , 3, 1, 1);
      book(_h_x     , 4, 1, 1);
      book(_h_xi    , 5, 1, 1);
      book(_wSum,"TMP/wSum");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t numParticles = cfs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      _wSum->fill();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // Thrusts
      MSG_DEBUG("Calculating thrust");
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      for (const Particle& p : cfs.particles()) {
        const Vector3 mom3  = p.p3();
        const double energy = p.E();
        const double pTinT  = dot(mom3, thrust.thrustMajorAxis());
        const double pToutT = dot(mom3, thrust.thrustMinorAxis());
      	_h_pT_in ->fill(fabs(pTinT/GeV) );
      	_h_pT_out->fill(fabs(pToutT/GeV));
        const double momT = dot(thrust.thrustAxis(), mom3);
        const double rapidityT = 0.5 * std::log((energy + momT) / (energy - momT));
      	_h_y->fill(fabs(rapidityT));
        const double mom = mom3.mod();
        const double scaledMom = mom/meanBeamMom;
        const double logInvScaledMom = -std::log(scaledMom);
        _h_xi->fill(logInvScaledMom);
        _h_x ->fill(scaledMom      );
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_pT_in , 1./ *_wSum);
      scale(_h_pT_out, 1./ *_wSum);
      scale(_h_y     , 1./ *_wSum);
      scale(_h_x     , 1./ *_wSum);
      scale(_h_xi    , 1./ *_wSum);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pT_in,_h_pT_out,_h_y,_h_x,_h_xi;
    CounterPtr _wSum;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2003_I595335);


}
