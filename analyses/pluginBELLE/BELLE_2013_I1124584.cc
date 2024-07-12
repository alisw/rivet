// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Bs0 -> Ds* Ds*
  class BELLE_2013_I1124584 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2013_I1124584);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==531);
      declare(ufs, "UFS");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles BS0 = apply<UnstableParticles>(event, "UFS").particles();
      for(const Particle & p : BS0) {
	if(p.children().size()!=2) continue;
	if(p.children()[0].pid()!=-p.children()[1].pid()) continue;
	if(p.children()[0].abspid()!=433) continue;
	Particle Dp = p.children()[0];
	Particle Dm = p.children()[1];
	if     (p.pid()>0 && Dp.pid()<0) swap(Dp,Dm);
	else if(p.pid()<0 && Dp.pid()>0) swap(Dp,Dm);
	// boost to rest frame
	LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pB =  boostB.transform(p.momentum());
	if(Dp.children().size()==2) {
	  Particle gamma;
	  bool found = true;
	  if(Dp.children()[0].pid()==PID::GAMMA &&
	     Dp.children()[1].abspid()==431)
	    gamma = Dp.children()[0];
	  else if (Dp.children()[1].pid()==PID::GAMMA &&
		   Dp.children()[0].abspid()==431)
	    gamma = Dp.children()[1];
	  else
	    found = false;
	  if( found) {
	    FourMomentum pD     = boostB.transform(Dp.momentum());
	    FourMomentum pgamma = boostB.transform(gamma.momentum());
	    LorentzTransform boostD = LorentzTransform::mkFrameTransformFromBeta(pD.betaVec());
	    Vector3 axisB = boostD.transform(pB).p3().unit();
	    Vector3 axisG = boostD.transform(pgamma).p3().unit();
	    _h[0]->fill(axisB.dot(axisG));
	  }
	}
	if(Dm.children().size()==2) {
	  Particle gamma;
	  bool found = true;
	  if(Dm.children()[0].pid()==PID::GAMMA &&
	     Dm.children()[1].abspid()==431)
	    gamma = Dm.children()[0];
	  else if (Dm.children()[1].pid()==PID::GAMMA &&
		   Dm.children()[0].abspid()==431)
	    gamma = Dm.children()[1];
	  else
	    found = false;
	  if( found) {
	    FourMomentum pD     = boostB.transform(Dm.momentum());
	    FourMomentum pgamma = boostB.transform(gamma.momentum());
	    LorentzTransform boostD = LorentzTransform::mkFrameTransformFromBeta(pD.betaVec());
	    Vector3 axisB = boostD.transform(pB).p3().unit();
	    Vector3 axisG = boostD.transform(pgamma).p3().unit();
	    _h[1]->fill(axisB.dot(axisG));
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2013_I1124584);

}
