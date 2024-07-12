// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> X_u ell nu
  class BABAR_2012_I1079912 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2012_I1079912);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==511 or
				Cuts::abspid==521),"UFS");
      // histograms
      book(_h,1,1,1);
      book(_nB,"/TMP/nB");
    }

    void findDecayProducts(Particle parent, Particles & em, Particles & ep,
			   Particles & nue, Particles & nueBar, bool & charm) {
      for(const Particle & p : parent.children()) {
	if(PID::isCharmHadron(p.pid())) {
	  charm=true;
	}
	else if(p.pid() == PID::EMINUS) {
	  em.push_back(p);
	}
	else if(p.pid() == PID::EPLUS) {
	  ep.push_back(p);
	}
	else if(p.pid() == PID::NU_E || p.pid()==PID::NU_MU) {
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
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle& p : ufs.particles()) {
	if(p.children().empty()) continue;
	if(p.children()[0].abspid()==p.abspid()) continue;
	_nB->fill();
	bool charm = false;
	Particles em,ep,nue,nueBar;
	findDecayProducts(p,em,ep,nue,nueBar,charm);
	if(charm) continue;
	FourMomentum pl,pnu;
	if(em.size()==1 && nueBar.size()==1) {
	  pl  = em[0].momentum();
	}
	else if(ep.size()==1 && nue.size()==1) {
	  pl  = ep[0].momentum();
	}
	else
	  continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	pl  = boost.transform(pl );
	double pp = pl.p3().mod();
	for(const auto & bin : _h->bins())
	  if(bin.xMin()<pp) _h->fill(bin.xMid());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1e3 due br normalization and 0.1 to remove bin width
      scale(_h, 0.1*1e3/ *_nB);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _nB;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2012_I1079912);

}
