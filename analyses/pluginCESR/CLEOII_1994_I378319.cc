// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_1+ and D_2 spectra
  class CLEOII_1994_I378319 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1994_I378319);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_D2_cTheta,2,1,1);
      book(_h_D2_x     ,4,1,1);
      book(_h_D1_cTheta,3,1,1);
      book(_h_D1_x     ,4,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idD1 = 10413;
      static const int idD2 = 415;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idD1 or Cuts::abspid==idD2)) {
	if(p.abspid()==idD1) {
	  // spectrum
	  double Pmax = sqrt(sqr(Emax)-sqr(2.421));
	  double xp = p.momentum().p3().mod()/Pmax;
	  _h_D1_x->fill(xp);
	}
	else {
	  double Pmax = sqrt(sqr(Emax)-sqr(2.461));
	  double xp = p.momentum().p3().mod()/Pmax;
	  _h_D2_x->fill(xp);
	}
	int sign = p.pid()/p.abspid();
	Particle Dstar;
	if(p.children().size()!=2) continue;
	if(p.children()[0].pid()==sign*423 &&
	   p.children()[1].pid()==sign*211) {
	  Dstar = p.children()[0];
	}
	else if(p.children()[1].pid()==sign*423 &&
		p.children()[0].pid()==sign*211) {
	  Dstar = p.children()[1];
	}
	else {
	  continue;
	}
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
	// first boost to the D_1,2 rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pDstar = boost1.transform(Dstar.momentum());
	FourMomentum pPion  = boost1.transform(pion .momentum());
	// to D* rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	Vector3 axis = pDstar.p3().unit();
	FourMomentum pp = boost2.transform(pPion);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(p.abspid()==idD1)
	  _h_D1_cTheta->fill(cTheta);
	else
	  _h_D2_cTheta->fill(cTheta);
      }      
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_D2_cTheta);
      normalize(_h_D2_x     );
      normalize(_h_D1_cTheta);
      normalize(_h_D1_x     );     
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_D2_cTheta,_h_D2_x,_h_D1_cTheta,_h_D1_x;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1994_I378319);

}
