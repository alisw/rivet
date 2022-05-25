// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Omega decay asymmetries
  class BABAR_2006_I719581 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I719581);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // Book histograms
      book(_h_ctheta_xic   ,1,1,1);
      book(_h_ctheta_omegac,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Xi_c0 baryons and Omega_c0 baryons
      for(const Particle& baryon : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4132 || Cuts::abspid==4332 )) {
	int sign = baryon.pid()/baryon.abspid();
	if(baryon.children().size()!=2) continue;
	Particle baryon1,meson1;
	if(baryon.abspid()==4132) {
	  if(baryon.children()[0].pid()==sign*3334 && 
	     baryon.children()[1].pid()==sign*321) {
	    baryon1 = baryon.children()[0];
	    meson1  = baryon.children()[1];
	  }
	  else if(baryon.children()[1].pid()==sign*3332 && 
		  baryon.children()[0].pid()==sign*321) {
	    baryon1 = baryon.children()[1];
	    meson1  = baryon.children()[0];
	  }
	  else
	    continue;
	}
	else {
	  if(baryon.children()[0].pid()==sign*3334 && 
	     baryon.children()[1].pid()==sign*211) {
	    baryon1 = baryon.children()[0];
	    meson1  = baryon.children()[1];
	  }
	  else if(baryon.children()[1].pid()==sign*3334 && 
		  baryon.children()[0].pid()==sign*211) {
	    baryon1 = baryon.children()[1];
	    meson1  = baryon.children()[0];
	  }
	  else
	    continue;
	}
	Particle baryon2,meson2;
	if(baryon1.children()[0].pid()== sign*3122 && 
	   baryon1.children()[1].pid()==-sign*321) {
	  baryon2 = baryon1.children()[0];
	  meson2  = baryon1.children()[1];
	}
	else if(baryon1.children()[1].pid()== sign*3122 && 
		baryon1.children()[0].pid()==-sign*321) {
	  baryon2 = baryon1.children()[1];
	  meson2  = baryon1.children()[0];
	}
	else
	  continue;
	// first boost to the Xic/Omegac rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(baryon.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	// to omega rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(baryon.abspid()==4132)
	  _h_ctheta_xic->fill(cTheta,1.);
	else 
	  _h_ctheta_omegac->fill(cTheta,1.);	
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_ctheta_xic);
      normalize(_h_ctheta_omegac);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ctheta_xic,_h_ctheta_omegac;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2006_I719581);


}
