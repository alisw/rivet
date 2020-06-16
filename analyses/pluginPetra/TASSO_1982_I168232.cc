// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi0 spectrum at 14 and 34 GeV
  class TASSO_1982_I168232 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1982_I168232);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {


      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      
      // Book histograms
      if(fuzzyEquals(sqrtS()/GeV, 14., 1e-3)) {
	book(_h_E, 2,1,1);
	book(_h_p, 2,2,2);
	book(_h_x, 2,3,3);
      }
      else if (fuzzyEquals(sqrtS()/GeV, 34., 1e-3)) {
	book(_h_E, 3,1,1);
	book(_h_p, 3,2,2);
	book(_h_x, 3,3,3);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::PI0)) {
	if(!p.parents().empty() && p.parents()[0].pid()==PID::K0S)
	  continue;
	double xE = p.E()/meanBeamMom;
	double beta = p.p3().mod()/p.E();
	_h_E->fill(p.E()       );
	_h_p->fill(p.p3().mod());
	_h_x->fill(xE,1./beta);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_E, crossSection()/nanobarn/sumOfWeights());
      scale(_h_p, crossSection()/nanobarn/sumOfWeights());
      scale(_h_x, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_E,_h_p,_h_x;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1982_I168232);

}
