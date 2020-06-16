// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Electron spectrum in B decays
  class BABAR_2017_I1498564 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2017_I1498564);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // Book histograms
      // specify custom binning
      book(_h_light, 1,1,1);
      book(_h_all  , 1,1,2);
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
	else if(p.pid() == PID::NU_E) {
	  nue.push_back(p);
	}
	else if(p.pid() == PID::NU_EBAR) {
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
      // find and loop over Upslion(4S)
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==300553)) {
      	// boost to rest frame
      	LorentzTransform cms_boost;
      	if (p.p3().mod() > 1*MeV)
      	  cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
      	// loop over decay products
      	for(const Particle & p : p.children()) {
      	  if(p.abspid()==511 || p.abspid()==521) {
	    _nB->fill();
	    bool charm = false;
      	    Particles em,ep,nue,nueBar;
	    findDecayProducts(p,em,ep,nue,nueBar,charm);
	    if(em.size()==1 && nueBar.size()==1) {
	      double pmod = cms_boost.transform(em[0].momentum()).p3().mod();
	      _h_all->fill(pmod);
	      if(!charm) _h_light->fill(pmod);
	    }
	    else if(ep.size()==1 && nue.size()==1) {
	      double pmod = cms_boost.transform(ep[0].momentum()).p3().mod();
	      _h_all->fill(pmod);
	      if(!charm) _h_light->fill(pmod);
	    }
      	  }
      	}	
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalize to number of B decays
      scale(_h_all, 1./ *_nB);
      scale(_h_light, 1./ *_nB);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_light,_h_all;
    CounterPtr _nB;
    //@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2017_I1498564);

}
