// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CELLO_1982_I12010 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1982_I12010);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {      
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      unsigned int iloc(0);
      if(isCompatibleWithSqrtS(22.)) {
	iloc=1;
      }
      else if (isCompatibleWithSqrtS(34.)) {
	iloc=2;
      }
      else
	MSG_ERROR("Beam energy not supported!");
      book(_histEEC  , 1, 1, iloc);
      book(_histAEEC , 3, 1, iloc);
      book(_weightSum, "TMP/weightSum");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if ( fs.particles().size() < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      _weightSum->fill();

      double Evis = 0.0;
      for (const Particle& p : fs.particles()) {
        Evis += p.E();
      }
      double Evis2 = sqr(Evis);
      // (A)EEC
      // Need iterators since second loop starts at current outer loop iterator, i.e. no "foreach" here!
      for (Particles::const_iterator p_i = fs.particles().begin(); p_i != fs.particles().end(); ++p_i) {
        for (Particles::const_iterator p_j = p_i; p_j != fs.particles().end(); ++p_j) {
          const Vector3 mom3_i = p_i->momentum().p3();
          const Vector3 mom3_j = p_j->momentum().p3();
          const double energy_i = p_i->momentum().E();
          const double energy_j = p_j->momentum().E();
          const double thetaij = mom3_i.unit().angle(mom3_j.unit());
          double eec = (energy_i*energy_j) / Evis2;
	  if(p_i != p_j) eec *= 2.;
          _histEEC->fill(thetaij, eec);
          if (thetaij <0.5*M_PI)
            _histAEEC->fill( thetaij, -eec);
          else
            _histAEEC->fill( M_PI-thetaij, eec);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histEEC,  1.0/ *_weightSum);
      scale(_histAEEC, 1.0/ *_weightSum);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _histEEC, _histAEEC;
    CounterPtr _weightSum;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CELLO_1982_I12010);


}
