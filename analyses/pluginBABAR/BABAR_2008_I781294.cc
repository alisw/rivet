// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi* decay asymmetries
  class BABAR_2008_I781294 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2008_I781294);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // book histos
      book(_h_ctheta,1,1,1);
      book(_h_P0    ,2,1,1);
      book(_h_P2    ,2,1,2);
      book(_wgtSum,"/TMP/WSUM");
    }

    void findChildren(const Particle &p, int & sign,
		      unsigned int & nprod, Particles & Xi, Particles &pi,Particles &K) {
      for(const Particle & child : p.children()) {
	if(child.pid()==sign*3312) {
	  Xi.push_back(child);
	  ++nprod;
	}
	else if(child.pid()==211) {
	  pi.push_back(child);
	  ++nprod;
	}
	else if(child.pid()==321) {
	  K.push_back(child);
	  ++nprod;
	}
	else if(child.children().empty()) {
	  ++nprod;
	}
	else {
	  findChildren(child,sign,nprod,Xi,pi,K);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_c
      for(const Particle& baryon : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4122)) {
	int sign = baryon.pid()/baryon.abspid();
	Particles Xi,pi,K;
	unsigned int nprod(0);
	findChildren(baryon,sign,nprod,Xi,pi,K);
	// check Lambda_c -> Xi- pi+ K+ decya mode
	if(nprod!=3||Xi.size()!=1||pi.size()!=1||K.size()!=1) continue;
	// first boost to the Lambda_c rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(baryon.momentum().betaVec());
	FourMomentum pbaryon1 = Xi[0].momentum()+pi[0].momentum();
	pbaryon1              = boost1.transform(pbaryon1);
	FourMomentum pbaryon2 = boost1.transform(Xi[0].momentum());
	// to Xi* frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	double mass = pbaryon1.mass();
	if(mass>1.5 && mass< 1.65)
	  _h_ctheta->fill(cTheta);
	_h_P0->fill(mass);
	_h_P2->fill(mass,0.5*(3.*sqr(cTheta)-1.));
	_wgtSum->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_ctheta);
      normalize(_h_P0,1.,false);
      normalize(_h_P2,1.,false);
    }
    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_ctheta, _h_P0,_h_P2;
    CounterPtr _wgtSum;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2008_I781294);

}
