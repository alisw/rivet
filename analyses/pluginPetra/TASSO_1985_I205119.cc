// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief K0 and Lambda spectra at 14, 22 and 34 GeV.
  class TASSO_1985_I205119 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1985_I205119);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      sqs = 1.;
      if(isCompatibleWithSqrtS(14.)) {
	book(_h_kaon_x  ,  1,1,1);
	book(_h_lambda_x,  4,1,1);
	book(_h_kaon_p  ,  7,1,1);
	book(_h_lambda_p, 10,1,1);
	sqs = 14.;
      }
      else if (isCompatibleWithSqrtS(22.)) {
	book(_h_kaon_x  ,  2,1,1);
	book(_h_lambda_x,  5,1,1);
	book(_h_kaon_p  ,  8,1,1);
	book(_h_lambda_p, 11,1,1);
	sqs = 22.;
      }
      else if (isCompatibleWithSqrtS(34.)) {
	book(_h_kaon_x  , 3,1,1);
	book(_h_lambda_x, 6,1,1);
	book(_h_kaon_p  , 9,1,1);
	book(_h_lambda_p,12,1,1);
	sqs = 34.;
      }
      else
	MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      for (const Particle& p : apply<UnstableParticles>(event, "UFS").
	       particles(Cuts::abspid==PID::LAMBDA or Cuts::pid==130 or Cuts::pid==310)) {
	double xE = p.E()/meanBeamMom;
	double modp = p.p3().mod();
	double beta = modp/p.E();
	if(abs(p.pid())==PID::LAMBDA) {
	  _h_lambda_x->fill(xE,1./beta);
	  _h_lambda_p->fill(modp,1.);
	}
	else {
	  _h_kaon_x->fill(xE,1./beta);	 
	  _h_kaon_p->fill(modp,1.);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_kaon_x  , sqr(sqs)*crossSection()/microbarn/sumOfWeights());
      scale(_h_lambda_x, sqr(sqs)*crossSection()/microbarn/sumOfWeights());
      scale(_h_kaon_p  , crossSection()/nanobarn/sumOfWeights());
      scale(_h_lambda_p, crossSection()/nanobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_kaon_x, _h_lambda_x, _h_kaon_p, _h_lambda_p;
    double sqs;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1985_I205119);


}
