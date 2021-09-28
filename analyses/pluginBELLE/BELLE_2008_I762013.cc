// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s1 decay angles
  class BELLE_2008_I762013 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2008_I762013);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_alpha,1,1,1);
      book(_h_beta ,2,1,1);
      book(_h_gamma,3,1,1);
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
      	Particle Dstar,Kaon;
      	if(p.children().size()!=2) continue;
      	if(p.children()[0].pid()==sign*413 &&
      	   isK0(p.children()[1].pid())) {
      	  Dstar = p.children()[0];
      	  Kaon = p.children()[1];
      	}
      	else if(p.children()[1].pid()==sign*413 &&
      		isK0(p.children()[0].pid())) {
      	  Kaon = p.children()[0];
      	  Dstar = p.children()[1];
      	}
      	else {
      	  continue;
      	}
      	// first boost to the D_s1 rest frame
      	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
      	FourMomentum pDstar = boost1.transform(Dstar.momentum());
      	double cTheta = pDstar.p3().unit().dot(p.momentum().p3().unit());
      	_h_alpha->fill(cTheta);
      	if(Dstar.children().size()!=2) continue;
      	Particle D0,Pion;
      	if(Dstar.children()[0].pid()== sign*211 && 
      	   Dstar.children()[1].pid()== sign*421) {
      	  Pion = Dstar.children()[0];
      	  D0   = Dstar.children()[1];
      	}
      	else if(Dstar.children()[1].pid()== sign*211 && 
      		Dstar.children()[0].pid()== sign*421) {
      	  D0   = Dstar.children()[0];
      	  Pion = Dstar.children()[1];
      	}
      	else
      	  continue;
      	// boost to D_s frame
      	FourMomentum pD  = boost1.transform(D0  .momentum());
      	FourMomentum pK  = boost1.transform(Kaon.momentum());
      	FourMomentum pPi = boost1.transform(Pion.momentum());
      	// to D* rest frame
      	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
      	Vector3 axis = pDstar.p3().unit();
      	FourMomentum pp = boost2.transform(pD);
      	// calculate angle
      	double cThetap = pp.p3().unit().dot(axis);
      	_h_gamma->fill(cThetap);
	// finally beta
	Vector3 n1 = pD  .p3().cross(pPi.p3()).unit();
	Vector3 n2 = axis.cross(pK.p3()).unit();
	double beta = acos(n1.dot(n2));
	_h_beta->fill(beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_alpha);
      normalize(_h_beta );
      normalize(_h_gamma);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_alpha, _h_beta, _h_gamma;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2008_I762013);

}
