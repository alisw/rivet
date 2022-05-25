// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief EEC at 29 GeV
  class MARKII_1988_I250899 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1988_I250899);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      for(unsigned int ix=0;ix<3;++ix) {
	book(_h_EEC [ix],2*ix+1,1,1);
	book(_h_AEEC[ix],2*ix+2,1,1);
      }
      book(_weightSum, "TMP/weightSum");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      if ( fs.particles().size() < 2) vetoEvent;
      _weightSum->fill();
      // visible energy
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
          const double thetaij = 180.*mom3_i.unit().angle(mom3_j.unit())/M_PI;
          double eec = (energy_i*energy_j) / Evis2;
	  if(p_i != p_j) eec *= 2.;
          for(unsigned int ix=0;ix<3;++ix) _h_EEC[ix]->fill(thetaij, eec);
          if (thetaij <90.)
            for(unsigned int ix=0;ix<3;++ix) _h_AEEC[ix]->fill( thetaij, -eec);
          else
            for(unsigned int ix=0;ix<3;++ix) _h_AEEC[ix]->fill( 180.-thetaij, eec);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	// factor due angle in degress not radians
	scale(_h_EEC [ix], 180.0/M_PI/ *_weightSum);
	scale(_h_AEEC[ix], 180.0/M_PI/ *_weightSum);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_EEC[3],_h_AEEC[3];
    CounterPtr _weightSum;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1988_I250899);

}
