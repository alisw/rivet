// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Spectrum for D_s1 
  class CLEOII_1993_I352823 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1993_I352823);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_x     ,3,1,1);
      book(_h_cTheta,4,1,1);
      book(_r[0],2,1,1);
      book(_r[1],2,1,2);
    }

    bool isK0(int id) {
      return id==310 || id==130 || abs(id)==311;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int DsID = 10433;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.535));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==DsID)) {
	// spectrum
	double xp = p.momentum().p3().mod()/Pmax;
        _h_x->fill(xp);
	// decay angle
	int sign = p.pid()/DsID;
	Particle Dstar;
	if(p.children().size()!=2) continue;
	if(p.children()[0].pid()==sign*423 &&
	   p.children()[1].pid()==sign*321) {
	  Dstar = p.children()[0];
	}
	else if(p.children()[1].pid()==sign*423 &&
		p.children()[0].pid()==sign*321) {
	  Dstar = p.children()[1];
	}
	else if(p.children()[0].pid()==sign*413 &&
		isK0(p.children()[1].pid())) {
	  _r[1]->fill(0.5);
	  continue;
      	}
      	else if(p.children()[1].pid()==sign*413 &&
      		isK0(p.children()[0].pid())) {
	  _r[1]->fill(0.5);
	  continue;
      	}
	else {
	  continue;
	}
	_r[0]->fill(0.5);
	if(Dstar.children().size()!=2) continue;
	Particle pion;
	if(Dstar.children()[0].pid()== 111 && 
	   Dstar.children()[1].pid()== sign*421) {
	  pion = Dstar.children()[0];
	}
	else if(Dstar.children()[1].pid()== 111 && 
		Dstar.children()[0].pid()== sign*421) {
	  pion = Dstar.children()[1];
	}
	else
	  continue;
	// first boost to the D_s1 rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pDstar = boost1.transform(Dstar.momentum());
	FourMomentum pPion  = boost1.transform(pion .momentum());
	// to D* rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	Vector3 axis = pDstar.p3().unit();
	FourMomentum pp = boost2.transform(pPion);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h_cTheta->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x     );
      normalize(_h_cTheta);
      scale(_r[0],crossSection()/sumOfWeights()/picobarn);
      scale(_r[1],crossSection()/sumOfWeights()/picobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x,_h_cTheta;
    Histo1DPtr _r[2];
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1993_I352823);

}
