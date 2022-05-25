// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief K+/K0 specta at 29 GeV
  class MARKII_1985_I207785 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1985_I207785);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      //Histograms
      book(_h_K0,2,1,1);
      book(_h_Kp,4,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==321 or Cuts::abspid==310 or Cuts::abspid==130)) {
	const double xp = 2.*p.p3().mod()/sqrtS();
	const double beta = p.p3().mod() / p.E();
	if(p.abspid()==321)
	  _h_Kp->fill(xp,1./beta);
	else
	  _h_K0->fill(xp,1./beta);
      }
    }
    

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_K0, 1./sumOfWeights());
      scale(_h_Kp, 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_K0,_h_Kp;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1985_I207785);

}
