// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0 and eta spectra at 91 GeV
  class L3_1994_I374698 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_1994_I374698);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_pi0_x_P , 5,1,1);
      book(_h_pi0_x_H , 6,1,1);
      book(_h_pi0_xi  , 7,1,1);
      book(_h_eta_x_P , 8,1,1);
      book(_h_eta_x_H , 9,1,1);
      book(_h_eta_xi  ,10,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle & p : ufs.particles(Cuts::pid==PID::PI0 ||
					      Cuts::pid==PID::ETA)) {
        const double x = 2.*p.momentum().p3().mod()/sqrtS();
	if(p.pid()==PID::PI0) {
	  _h_pi0_x_P->fill(x);
	  _h_pi0_x_H->fill(x);
	  _h_pi0_xi->fill(-log(x));
	}
	else {
	  _h_eta_x_P->fill(x);
	  _h_eta_x_H->fill(x);
	  _h_eta_xi->fill(-log(x));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pi0_x_P   , 1./sumOfWeights());
      scale(_h_pi0_x_H   , 1./sumOfWeights());
      scale(_h_pi0_xi    , 1./sumOfWeights());
      scale(_h_eta_x_P   , 1./sumOfWeights());
      scale(_h_eta_x_H   , 1./sumOfWeights());
      scale(_h_eta_xi    , 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pi0_x_P,_h_pi0_x_H,_h_pi0_xi;
    Histo1DPtr _h_eta_x_P,_h_eta_x_H,_h_eta_xi;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(L3_1994_I374698);

}
