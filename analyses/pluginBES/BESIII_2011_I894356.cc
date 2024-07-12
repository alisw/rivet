// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief chi_c1 -> gamma +V
  class BESIII_2011_I894356 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2011_I894356);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==20443), "UFS");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
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
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles()) {
	if(p.children().size()!=2) continue;
	Particle gamma,meson;
	// get the photon and the meson
	if(p.children()[0].pid()==PID::PHOTON) {
	  gamma = p.children()[0];
	  meson = p.children()[1];
	}
	else if(p.children()[1].pid()==PID::PHOTON) {
	  gamma = p.children()[1];
	  meson = p.children()[0];
	}
	else
	  continue;
	// check the meson
	int imode=-1;
	if(meson.pid()==333)      imode = 0;
	else if(meson.pid()==113) imode = 1;
	else if(meson.pid()==223) imode = 2;
	else continue;
	Particle d1,d2,d3;
	if(imode==0) {
	  if(meson.children().size()!=2) continue;
	  if(meson.children()[0].pid()==-meson.children()[1].pid() &&
	     meson.children()[0].abspid()==321) {
	    d1 = meson.children()[0];
	    d2 = meson.children()[0];
	  }
	  else
	    continue;
	}
	else if(imode==1) {
	  if(meson.children().size()!=2) continue;
	  if(meson.children()[0].pid()==-meson.children()[1].pid() &&
	     meson.children()[0].abspid()==211) {
	    d1 = meson.children()[0];
	    d2 = meson.children()[0];
	  }
	  else
	    continue;
	}
	else if(imode==2) {
	  unsigned int ncount=0;
	  Particles pip,pim,pi0;
	  findChildren(meson,pim,pip,pi0,ncount);
	  if(ncount==3 && pim.size()==1 && pip.size()==1 && pi0.size()==1) {
	    d1=pip[0];
	    d2=pim[0];
	    d3=pi0[0];
	  }
	  else continue;
	}
	if(d1.pid()<0) swap(d1,d2);
	// boost in chi_c1 frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pMeson = boost1.transform(meson.momentum());
	FourMomentum p1=boost1.transform(d1.momentum());
	FourMomentum p2=boost1.transform(d2.momentum());
	Vector3 axis1 = pMeson.p3().unit();
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pMeson.betaVec());
	p1=boost2.transform(p1);
	p2=boost2.transform(p2);
	// axis in meson rest frame
	Vector3 axis2;
	if(imode<2) {
	  axis2=p1.p3().unit();
	}
	else {
	  FourMomentum p3 = boost1.transform(d3.momentum());
	  p3 = boost2.transform(p3);
	  axis2=p1.p3().cross(p2.p3()).unit();
	}
	_h[imode]->fill(axis1.dot(axis2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2011_I894356);

}
