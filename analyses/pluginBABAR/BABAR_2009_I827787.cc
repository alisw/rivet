// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> Xc ell - nu_ell
  class BABAR_2009_I827787 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I827787);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==511 ||
				Cuts::abspid==521),"UFS");
      // histograms
      for(unsigned int ix=0;ix<6;++ix) {
	book(_p_X[ix],1,1,1+ix);
	if(ix<3) book(_p_n[ix],2,1,1+ix);
      }
    }
    void findDecayProducts(Particle parent, Particles & em, Particles & ep,
			   Particles & nue, Particles & nueBar, bool & charm) {
      for(const Particle & p : parent.children()) {
	if(PID::isCharmHadron(p.pid())) {
	  charm=true;
	}
	else if(p.pid() == PID::EMINUS || p.pid()==PID::MUON) {
	  em.push_back(p);
	}
	else if(p.pid() == PID::EPLUS || p.pid()==PID::ANTIMUON) {
	  ep.push_back(p);
	}
	else if(p.pid() == PID::NU_E  || p.pid()==PID::NU_MU) {
	  nue.push_back(p);
	}
	else if(p.pid() == PID::NU_EBAR || p.pid()==PID::NU_MUBAR) {
	  nueBar.push_back(p);
	}
	else if(PID::isBottomHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar,charm);
	}
	else if(!PID::isHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar,charm);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const double Lambda=0.65*GeV;
      // find and loop over Upslion(4S)
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	if(p.children().empty() ||
	   (p.children().size()==1 && p.children()[1].abspid()==p.abspid()))
	  continue;
	// find decay products
	bool charm = false;
	Particles em,ep,nue,nueBar;
	findDecayProducts(p,em,ep,nue,nueBar,charm);
	if(!charm) continue;
	FourMomentum pl,pnu;
	if(em.size()==1 && nueBar.size()==1 && em[0].pid()+1==-nueBar[0].pid()) {
	  pl  = em[0].momentum();
	  pnu = nueBar[0].momentum();
	}
	else if(ep.size()==1 && nue.size()==1 && nue[0].pid()==-ep[0].pid()+1) {
	  pl  = ep[0].momentum();
	  pnu = nue[0].momentum();
       	}
	else
	  continue;
	// boost to rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pX = boost.transform(p.momentum()-pl-pnu);
	pl = boost.transform(pl);
	double modp = pl.p3().mod();
	double mX  = pX.mass();
	double nX2 = sqr(mX)-2.*Lambda*pX.t()+sqr(Lambda);
	double m1=mX,m2=nX2;
	for(unsigned int ix=0;ix<6;++ix) {
	  for(const auto & bin : _p_X[ix]->bins())
	    if(bin.xMin()<modp) _p_X[ix]->fill(bin.xMid(),m1);
	  m1 *=mX;
	  if(ix>=3) continue;
	  for(const auto & bin : _p_n[ix]->bins())
	    if(bin.xMin()<modp) _p_n[ix]->fill(bin.xMid(),m2);
	  m2 *=nX2;
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
    }

    /// @}


    /// @name Histograms
    /// @{
    Profile1DPtr _p_X[6],_p_n[3];

    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2009_I827787);

}
