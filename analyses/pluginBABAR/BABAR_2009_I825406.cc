// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B0 -> a1+ a1-
  class BABAR_2009_I825406 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I825406);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==511);
      declare(ufs, "B0");
      // histograms
      book(_p,1,1,1);
    }

    void findChildren(const Particle & p, Particles & pim, Particles & pip,
		      Particles & pi0, unsigned int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  pip.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PIMINUS) {
	  pim.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PI0) {
	  pi0.push_back(child);
	  ncount+=1;
	}
	else if(child.children().empty()) {
	  ncount+=1;
	}
    	else
    	  findChildren(child,pim,pip,pi0,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles B0 = apply<UnstableParticles>(event, "B0").particles();
      for(const Particle & p : B0) {
	// skip cases with mixing
	if(p.children().size()==1 && p.children()[0].abspid()==p.abspid()) continue;
	// particle antiparticle pair of a_1
	if(p.children().size()!=2 || p.children()[0].pid()!=-p.children()[1].pid() ||
	   p.children()[0].abspid()!=20213) continue;
	Particle a1p = p.children()[0], a1m = p.children()[1];
	if( (p.pid()>0&&a1p.pid()<0) || (p.pid()<0&&a1p.pid()>0) ) swap(a1p,a1m);
	Particles pip,pim,pi0;
	unsigned int ncount=0;
	findChildren(a1p,pim,pip,pi0,ncount);
	if(ncount!=3 || pip.size()!=2 || pim.size()!=1) continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pa1p = boost.transform(a1p.momentum());
	FourMomentum pa1m = boost.transform(a1m.momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pa1p.betaVec());
	FourMomentum ppip = boost2.transform(boost.transform(pip[0].momentum()));
	FourMomentum ppim = boost2.transform(boost.transform(pim[0].momentum()));
	Vector3 n = ppip.p3().cross(ppim.p3()).unit();
	double cTheta = n.dot(pa1m.p3().unit());
	_p->fill(5.28,(2.-5.*sqr(cTheta)));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I825406);

}
