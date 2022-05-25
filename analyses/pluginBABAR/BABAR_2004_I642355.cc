// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief M_x in B -> eta' Xs
  class BABAR_2004_I642355 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2004_I642355);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histograms
      book(_w_ups4,"/TMP/w_ups4");
      book(_h_X,1,1,1);
    }

    void findDecay(Particle parent, bool & charm, Particles & etaP, Particles & pi0,
		   Particles & pip, Particles & Kp, Particles & K0) {
      for(const Particle & p : parent.children()) {
	if(isCharm(p)) {
	  charm = true;
	  return;
	}
	else if(p.pid()==PID::ETAPRIME) {
	  etaP.push_back(p);
	}
	else if(p.abspid()==PID::PIPLUS) {
	  pip.push_back(p);
	}
	else if(p.abspid()==PID::KPLUS) {
	  Kp.push_back(p);
	}
	else if(p.pid()==PID::PI0) {
	  pi0.push_back(p);
	}
	else if(p.pid()==PID::K0S) {
	  K0.push_back(p);
	}
	else if(!p.children().empty()) {
	  findDecay(p,charm,etaP,pi0,pip,Kp,K0); 
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle & p : ufs.particles(Cuts::pid==300553)) {
	_w_ups4->fill();
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
	for(const Particle & B:p.children()) {
	  if(B.abspid()!=511 && B.abspid()!=521) continue;
	  // find the decay
	  bool charm=false;
	  Particles etaP, pi0, pip, Kp, K0;
	  findDecay(B,charm,etaP,pi0,pip,Kp,K0);
	  if(charm || etaP.size()!=1 ) continue;
	  if(Kp.empty()&&Kp.empty()) continue;
	  double pEta = boost.transform(etaP[0].momentum()).p3().mod();
	  if(pEta<2. || pEta>2.7) continue;
	  // four momentum of the Xs system
	  FourMomentum pX;
	  if(K0.size()==1)
	    pX+= K0[0].momentum();
	  else if(Kp.size()==1)
	    pX+= Kp[0].momentum();
	  else
	    continue;
	  for(const Particle & pi : pip)
	    pX+=pi.momentum();
	  if(pi0.size()==1)
	    pX+=pi0[0].momentum();
	  _h_X->fill(pX.mass());
	  
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_X, 0.5/ *_w_ups4);
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _w_ups4;
    Histo1DPtr _h_X;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2004_I642355);

}
