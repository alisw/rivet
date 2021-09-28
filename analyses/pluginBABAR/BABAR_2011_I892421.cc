// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s1 decay angles
  class BABAR_2011_I892421 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2011_I892421);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_cP,1,1,1);
      book(_h_c ,2,1,1);
    }

    bool isK0(int id) {
      return id==310 || id==130 || abs(id)==311;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int DsID = 10433;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==DsID)) {
	// decay angle
	int sign = p.pid()/DsID;
	Particle Dstar;
	if(p.children().size()!=2) continue;
	if(p.children()[0].pid()==sign*413 &&
	   isK0(p.children()[1].pid())) {
	  Dstar = p.children()[0];
	}
	else if(p.children()[1].pid()==sign*413 &&
		isK0(p.children()[0].pid())) {
	  Dstar = p.children()[1];
	}
	else {
	  continue;
	}
	// first boost to the D_s1 rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pDstar = boost1.transform(Dstar.momentum());
	double cTheta = pDstar.p3().unit().dot(p.momentum().p3().unit());
	_h_c->fill(cTheta);
	if(Dstar.children().size()!=2) continue;
	Particle D0;
	if(Dstar.children()[0].pid()== sign*211 && 
	   Dstar.children()[1].pid()== sign*421) {
	  D0 = Dstar.children()[1];
	}
	else if(Dstar.children()[1].pid()== sign*211 && 
		Dstar.children()[0].pid()== sign*421) {
	  D0 = Dstar.children()[0];
	}
	else
	  continue;
	// boost to D_s frame
	FourMomentum pD  = boost1.transform(D0.momentum());
	// to D* rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	Vector3 axis = pDstar.p3().unit();
	FourMomentum pp = boost2.transform(pD);
	// calculate angle
	double cThetap = pp.p3().unit().dot(axis);
	_h_cP->fill(cThetap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_cP);
      normalize(_h_c );
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_cP,_h_c;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2011_I892421);

}
