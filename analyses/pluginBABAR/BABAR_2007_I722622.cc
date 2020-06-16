// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c' spectrum
  class BABAR_2007_I722622 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2007_I722622);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableFinalState(),"UFS");
      // Book histograms
      book(_h_p_0,3,1,2);
      book(_h_p_p,3,1,1);
      book(_h_ctheta,4,1,1);
    }

    void findChildren(Particle parent, unsigned int & nStable,
		      Particles &Xi, unsigned int &nPi) {
      for(const Particle & p : parent.children()) {
	if(p.abspid()==PID::PIPLUS) {
	  ++nPi;
	  ++nStable;
	}
	else if(p.abspid()==PID::XIMINUS) {
	  Xi.push_back(p);
	  ++nStable;
	}
	else if(!p.children().empty()) {
	  findChildren(p,nStable,Xi,nPi);
	}
	else {
	  ++nStable;
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int id0 = 4312, idp = 4322;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idp or Cuts::abspid==id0)) {
	int idXic=4232;
	if(p.abspid()==idp)
	  _h_p_p->fill(p.momentum().p3().mod());
	else {
	  _h_p_0->fill(p.momentum().p3().mod());
	  idXic=4132;
	}
	// first find the Xi_c in the decay
	if(p.children().size()!=2) continue;
	int sign=p.pid()/p.abspid();
	Particle Xi_c;
	if(p.children()[0].pid()==sign*idXic &&
	   p.children()[1].pid()==22) {
	  Xi_c = p.children()[0];
	}
	else if(p.children()[1].pid()==sign*idXic &&
		p.children()[0].pid()==22) {
	  Xi_c = p.children()[1];
	}
	else
	  continue;
	// and the children of the Xi_c
	Particles Xi;
	unsigned int nStable(0),nPi(0);
	findChildren(Xi_c,nStable,Xi,nPi);
	if(Xi.size()!=1) continue;
	if( !(p.abspid()==idp && nPi==2 && nStable==3) &&
	    !(p.abspid()==id0 && nPi==1 && nStable==2) )
	  continue;
	// boost to Xi'_c rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pXic = boost1.transform(Xi_c.momentum());
	FourMomentum pXi  = boost1.transform(Xi[0].momentum());
	// then to Xi_c rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pXic.betaVec());
	Vector3 axis = pXic.p3().unit();
	FourMomentum pp = boost2.transform(pXi);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h_ctheta->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_p_0   );
      normalize(_h_p_p   );
      normalize(_h_ctheta);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_p_0,_h_p_p,_h_ctheta;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2007_I722622);

}
