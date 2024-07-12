// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c -> xi pi pi
  class BELLE_2018_I1698390 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1698390);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(),"UFS");
      // histo
      book(_h,1,1,1);
    }

    void findDecay(int & sign, const Particle & parent, Particles & xi, Particles & pi,
		   unsigned int & nstable) {
      for(const Particle & p : parent.children()) {
	if(p.pid()==sign*3312) {
	  ++nstable;
	  xi.push_back(p);
	}
	else if(p.pid()==sign*211) {
	  ++nstable;
	  pi.push_back(p);
	}
	else if(p.children().empty())
	  ++nstable;
	else
	  findDecay(sign,p,xi,pi,nstable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>(event,"UFS").particles(Cuts::abspid==4232)) {
	int sign = p.pid()/p.abspid();
	Particles xi,pi;
	unsigned int nstable(0);
	findDecay(sign,p,xi,pi,nstable);
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	if(nstable==3&&xi.size()==1&&pi.size()==2) {
	  double p1 = boost.transform(pi[0].momentum()).p3().mod();
	  double p2 = boost.transform(pi[1].momentum()).p3().mod();
	  double mass = p1>p2 ?
	    (pi[1].momentum()+xi[0].momentum()).mass() :
	    (pi[0].momentum()+xi[0].momentum()).mass();
	  _h->fill(mass);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1698390);

}
