// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c(2970) decays
  class BELLE_2021_I1809180 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1809180);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // set the PDG code
      _pid = getOption<double>("PID", 204232);
      // projections
      declare(UnstableParticles(Cuts::abspid==_pid), "UFS");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
	// xp cut
	double xp = p.momentum().p3().mod()/sqrt(0.25*sqrtS()-p.mass2());
	if (xp<.7) continue;
	// find the decay products
	// first decay
	if(p.children().size()!=2) continue;
	Particle XiStar,pi1;
	if(p.children()[0].abspid()==4314 &&
	   p.children()[1].abspid()== 211) {
	  XiStar = p.children()[0];
	  pi1    = p.children()[1];
	}
	else if(p.children()[1].abspid()==4314 &&
		p.children()[0].abspid()== 211) {
	  XiStar = p.children()[1];
	  pi1    = p.children()[0];
	}
	else
	  continue;
	Vector3 axis1 = p.momentum().p3().unit();
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum ppi1 = boost.transform(pi1.momentum());
	_h[0]->fill(axis1.dot(ppi1.p3().unit()));
	Particle pi2;
	// second decay
	if(XiStar.children()[0].abspid()==4232 &&
	   XiStar.children()[1].abspid()== 211) {
	  pi2    = XiStar.children()[1];
	}
	else if(XiStar.children()[1].abspid()== 4232 &&
		XiStar.children()[0].abspid()== 211) {
	  pi2    = XiStar.children()[0];
	}
	else
	  continue;
	FourMomentum pXiStar = boost.transform(XiStar.momentum());
	FourMomentum ppi2    = boost.transform(pi2.momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pXiStar.betaVec());
	ThreeVector axis2 = pXiStar.p3().unit();
	ppi2 =  boost2.transform(ppi2);
	_h[1]->fill(axis2.dot(ppi2.p3().unit()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    int _pid;
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1809180);

}
