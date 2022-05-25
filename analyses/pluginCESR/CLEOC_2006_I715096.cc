// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Electron spectrum in D decays
  class CLEOC_2006_I715096 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2006_I715096);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // Book histograms
      // specify custom binning
      book(_h_Dp, 1,1,1);
      book(_h_D0, 1,1,2);
    }

    void findDecayProducts(Particle parent, Particles & em, Particles & ep,
			   Particles & nue, Particles & nueBar) {
      for(const Particle & p : parent.children()) {
	if(p.pid() == PID::EMINUS) {
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
	else if(PID::isCharmHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar);
	}
	else if(!PID::isHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find and loop over psi(3770)
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::pid==30443)) {
	// boost to rest frame
	LorentzTransform boost;
	if (p.p3().mod() > 1*MeV)
	  boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	// loop over decay products
	for(const Particle & p2 : p.children()) {
	  if(p2.abspid()==411 || p2.abspid()==421) {
	    Particles em,ep,nue,nueBar;
	    findDecayProducts(p2,em,ep,nue,nueBar);
	    if(em.size()==1 && nueBar.size()==1) {
	      double pmod = boost.transform(em[0].momentum()).p3().mod();
	      if(p2.abspid()==411) _h_Dp->fill(pmod);
	      else                 _h_D0->fill(pmod);
	    }
	    else if(ep.size()==1 && nue.size()==1) {
	      double pmod = boost.transform(ep[0].momentum()).p3().mod();
	      if(p2.abspid()==411) _h_Dp->fill(pmod);
	      else                 _h_D0->fill(pmod);
	    }
	  }	
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_D0);
      normalize(_h_Dp);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Dp,_h_D0;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2006_I715096);

}
