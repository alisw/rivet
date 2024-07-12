// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief rho0, K*=, K*0 spectra at 29 GeV
  class HRS_1989_I276948 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1989_I276948);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_x_rho   ,1,1,1);
      book(_h_x_Kstar0,2,1,1);
      book(_h_x_Kstarp,2,1,2);
      book(_h_sig_rho    ,3,1,1);
      book(_h_sig_Kstar0 ,4,1,1);
      book(_h_mult_rho   ,3,1,3);
      book(_h_mult_Kstar0,4,1,3);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for(const Particle & p : ufs.particles(Cuts::abspid==113 ||
					     Cuts::abspid==323 ||
					     Cuts::abspid==313)) {
      	Vector3 mom3 = p.p3();
        const double energy = p.E();
      	double modp = mom3.mod();
      	double beta = modp/energy;
      	double xE = 2.*modp/sqrtS();
	if(p.pid()==113) {
	  _h_x_rho->fill(xE,1./beta);
	  _h_sig_rho->fill(sqrtS());
	  _h_mult_rho->fill(sqrtS());
	}
	else if(p.pid()==313) {
	  _h_x_Kstar0->fill(xE,1./beta);
	  _h_sig_Kstar0->fill(sqrtS());
	  _h_mult_Kstar0->fill(sqrtS());
	}
	else {
	  _h_x_Kstarp->fill(xE,1./beta);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale( _h_x_rho   , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale( _h_x_Kstar0, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale( _h_x_Kstarp, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale( _h_sig_rho   , crossSection()/picobarn/sumOfWeights());
      scale( _h_sig_Kstar0, crossSection()/picobarn/sumOfWeights());
      scale( _h_mult_rho   , 1./sumOfWeights());
      scale( _h_mult_Kstar0, 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x_rho,_h_x_Kstar0,_h_x_Kstarp;
    Histo1DPtr _h_sig_rho,_h_sig_Kstar0;
    Histo1DPtr _h_mult_rho,_h_mult_Kstar0;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1989_I276948);

}
