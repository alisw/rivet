// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief EEC for a wide range of energies
  class PLUTO_1981_I156315 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PLUTO_1981_I156315);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      // Book histograms
      unsigned int iloc(0);
      if(fuzzyEquals(sqrtS()/GeV, 7.7, 1e-3)) {
	iloc=1;
      }
      else if(fuzzyEquals(sqrtS()/GeV, 9.4, 1e-3)) {
	iloc=2;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 12., 1e-3)) {
	iloc=3;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 13., 1e-3)) {
	iloc=4;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 17., 1e-3)) {
	iloc=5;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 22., 1e-3)) {
	iloc=6;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 27.6, 1e-3)) {
	iloc=7;
      }
      else if (inRange(sqrtS()/GeV,30,31.6)) {
	iloc=8;
      }
      else
	MSG_ERROR("Beam energy not supported!");
      // Book histograms
      book(_h_EEC, 1, 1, iloc);
      if(iloc==7||iloc==8) {
	book(_h_AEEC, 5, 1, 1);
	// _h_opposite = bookHisto1D(2, 1, 1);
      }
      else if(iloc==21 ||iloc==2)
	book(_h_AEEC,4, 1, 1);
      book(_weightSum,"TMP/weightSum");
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
          const double thetaij = mom3_i.unit().angle(mom3_j.unit())/M_PI*180.;
          double eec = (energy_i*energy_j) / Evis2;
	  if(p_i != p_j) eec *= 2.;
	  _h_EEC ->fill(thetaij,  eec);
	  // if(_h_opposite) _h_opposite ->fill(mom3_i.unit().dot(mom3_j.unit()),  eec);
	  if(_h_AEEC) {
	    if (thetaij < 90.) {
	      _h_AEEC->fill(thetaij, -eec);
	    }
	    else {
	      _h_AEEC  ->fill(180.-thetaij, eec);
	    }
	  }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_EEC , 360.0/M_PI/ *_weightSum);
      scale(_h_AEEC, 360.0/M_PI/ *_weightSum);
      // scale(_h_opposite, 2./ *_weightSum);

    }

    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_EEC, _h_AEEC, _h_opposite;
    CounterPtr _weightSum;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PLUTO_1981_I156315);


}
