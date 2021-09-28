// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0 and charged particle spectra
  class L3_1991_I314407 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_1991_I314407);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_pi0_x     ,3,1,1);
      book(_h_pi0_xi    ,4,1,1);
      book(_h_charged_x ,5,1,1);
      book(_h_charged_xi,6,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      for (const Particle& p : fs.particles()) {
        const double x = 2.*p.momentum().p3().mod()/sqrtS();
	_h_charged_x ->fill(x);
	_h_charged_xi->fill(-log(x));


      }
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle & p : ufs.particles(Cuts::pid==PID::PI0)) {
        const double x = 2.*p.momentum().p3().mod()/sqrtS();
	_h_pi0_x ->fill(x);
	_h_pi0_xi->fill(-log(x));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pi0_x     , 1./sumOfWeights());
      scale(_h_pi0_xi    , 1./sumOfWeights());
      scale(_h_charged_x , 1./sumOfWeights());
      scale(_h_charged_xi, 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pi0_x,_h_pi0_xi,_h_charged_x,_h_charged_xi;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(L3_1991_I314407);

}
