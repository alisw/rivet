// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Hadron Spectra in $e^+e^-$ collisions at 29 GeV
  class HRS_1987_I215848 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1987_I215848);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_pi ,2,1,1);
      book(_h_Kp ,3,1,1);
      book(_h_p  ,4,1,1);
      book(_h_K0 ,5,1,1);
      book(_h_lam,6,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & p : ufs.particles()) {
	double xE = 2.*p.E()/sqrtS();
	const double beta = p.p3().mod() / p.E();
	if(p.pid()==130 || p.pid()==310)
	  _h_K0  ->fill(xE,1./beta);
	else if(p.abspid()==321)
	  _h_Kp  ->fill(xE,1./beta);
	else if(p.abspid()==211)
	  _h_pi ->fill(xE,1./beta);
	else if(p.abspid()==2212)
	  _h_p  ->fill(xE,1./beta);
	else if(p.abspid()==3122)
	  _h_lam->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pi , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      scale(_h_Kp , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      scale(_h_p  , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      scale(_h_K0 , crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
      scale(_h_lam, crossSection()*sqr(sqrtS())/microbarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pi,_h_Kp,_h_p,_h_K0,_h_lam;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(HRS_1987_I215848);

}
