// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief EEC at 34.6 GeV
  class PLUTO_1985_I215869 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PLUTO_1985_I215869);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // Initialise and register projections
      declare(FinalState(), "FS");
      // Book histograms
      book(_h_EEC, 1, 1, 1);
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
          const double cosij = mom3_i.unit().dot(mom3_j.unit());
          //const double thetaij = mom3_i.unit().angle(mom3_j.unit())/M_PI*180.;
          double eec = (energy_i*energy_j) / Evis2;
	  if(p_i != p_j) eec *= 2.;
	  _h_EEC ->fill(cosij,  eec);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_EEC , 1./ *_weightSum);
    }
    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_EEC;
    CounterPtr _weightSum;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PLUTO_1985_I215869);


}
