// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Excited Lambda_c spectra
  class CLEOII_1994_I381696 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1994_I381696);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_2595,5,1,1);
      book(_h_2625,6,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int id2595 = 14122;
      static const int id2625 = 4124;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==id2595 or
					     Cuts::abspid==id2625)) {
	// spectrum
	if(p.abspid()==id2595) {
	  double Pmax = sqrt(sqr(Emax)-sqr(2.595));
	  double xp = p.momentum().p3().mod()/Pmax;
	  _h_2595->fill(xp);
	}
	else {
	  double Pmax = sqrt(sqr(Emax)-sqr(2.625));
	  double xp = p.momentum().p3().mod()/Pmax;
	  _h_2625->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_2595);
      normalize(_h_2625);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_2595,_h_2625;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1994_I381696);

}
