// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi0 spectrum at 14 and 34 GeV
  class TASSO_1982_I168232 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1982_I168232);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {


      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      
      // Book histograms
      sqs = 1.0;
      if(isCompatibleWithSqrtS(14.)) {
	book(_h_E, 2,1,1);
	book(_h_p, 2,2,2);
	book(_h_x, 2,3,3);
	sqs = 14.0;
      }
      else if (isCompatibleWithSqrtS(34.)) {
	book(_h_E, 3,1,1);
	book(_h_p, 3,2,2);
	book(_h_x, 3,3,3);
	sqs = 34.0;
      }
      else
        MSG_ERROR("Not compatible with energy " << sqrtS() << "GeV.");
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
      scale(_h_x, sqr(sqs)*crossSection()/microbarn/sumOfWeights());

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_E,_h_p,_h_x;
    double sqs;
    //@}

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1982_I168232);

}
