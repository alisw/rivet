// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> D*0 omega
  class BABAR_2011_I920989 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2011_I920989);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      // histos
      for(unsigned int ix=0;ix<4;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1,1+ix,1+iy);
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
      for(const Particle & B0 : apply<UnstableParticles>(event,"UFS").particles()) {
	if(B0.children().size()!=2) continue;
	Particle Dstar,omega;
	if(B0.children()[0].abspid()==423 &&
	   B0.children()[1].pid()==223) {
	  Dstar = B0.children()[0];
	  omega = B0.children()[1];
	}
	else if (B0.children()[1].abspid()==423 &&
		 B0.children()[0].pid()==223) {
	  Dstar = B0.children()[1];
	  omega = B0.children()[0];
	}
	else
	  continue;
	// check the no of decay products
	if(Dstar.children().size()!=2 || omega.children().size()!=3)
	  continue;
	// find the children of the D* meson
	Particle D0;
	if(Dstar.children()[0].pid()==111 &&
	   Dstar.children()[1].abspid()==421)
	  D0 = Dstar.children()[1];
	else if(Dstar.children()[1].pid()==111 &&
		Dstar.children()[0].abspid()==421)
	  D0 = Dstar.children()[0];
	else
	  continue;
	// children of the omega
	unsigned int ncount=0;
	Particles pip,pim,pi0;
	findChildren(omega,pim,pip,pi0,ncount);
	if( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) continue;
	// boost to B rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.momentum().betaVec());
	FourMomentum pDstar = boost1.transform(Dstar.momentum());
	FourMomentum pD0    = boost1.transform(D0   .momentum());
	FourMomentum pomega = boost1.transform(omega.momentum());
	FourMomentum pPip   = boost1.transform(pip[0].momentum());
	FourMomentum pPim   = boost1.transform(pim[0].momentum());
	// boost to D* frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	pD0 = boost2.transform(pD0);
	double c1 = pD0.p3().unit().dot(pDstar.p3().unit());
	// boost to omega frame
	LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pomega.betaVec());
	pPip   = boost3.transform(pPip);
	pPim   = boost3.transform(pPim);
	Vector3 axisOmega = pPip.p3().cross(pPim.p3()).unit();
	double c2 = pomega.p3().unit().dot(axisOmega);
	for(unsigned int ix=0;ix<4;++ix) {
	  _h[ix][0]->fill(c1);
	  _h[ix][1]->fill(c2);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2011_I920989);

}
