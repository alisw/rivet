// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Upsilon -> Upsilon pi+pi-
  class BELLE_2017_I1610301 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1610301);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      for(unsigned int ix=0;ix<4;++ix) {
	book(_h_mpipi [ix],1,1,ix+1);
	book(_h_ctheta[ix],2,1,ix+1);
      }
    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles& pi0, Particles & onium) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
       	else if (id == PID::PIPLUS) {
       	  pip.push_back(p);
       	  ++nstable;
       	}
       	else if (id == PID::PI0) {
       	  pi0.push_back(p);
       	  ++nstable;
       	}
	else if (abs(id)%1000==443 || abs(id)%1000==553) {
	  onium.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,pi0,onium);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==100553 or
										 Cuts::pid==200553 or
										 Cuts::pid==300553)) {
	unsigned int nstable(0);
	Particles pip, pim, pi0, onium;
	findDecayProducts(ups,nstable,pip,pim,pi0,onium);
	// pions and no of products
	if(pip.size()!=1 || pim.size() !=1 || nstable !=3) continue;
	// check for onium
	if(onium.size() !=1) continue;
	if(!(onium[0].pid()==553 || (ups.pid()==300553 && onium[0].pid()==100553))) continue;
	FourMomentum q = pip[0].momentum()+pim[0].momentum();
	// boost particles to upsilon rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	FourMomentum ppi  = boost1.transform(pip[0].momentum());
	FourMomentum pUps = boost1.transform(onium[0].momentum());
	q = boost1.transform(q);
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(q.betaVec());
	ppi = boost2.transform(ppi);
	pUps = boost2.transform(pUps);
	double ctheta = ppi.p3().unit().dot(pUps.p3().unit());
	if(ups.pid()==100553) {
	  _h_mpipi[0]->fill(q.mass());
	  _h_ctheta[0]->fill(ctheta);
	}
	else if(ups.pid()==200553) {
	  _h_mpipi[1]->fill(q.mass());
	  _h_ctheta[1]->fill(ctheta);
	}
	else if(ups.pid()==300553) {
	  if(onium[0].pid()==553) {
	    _h_mpipi[2]->fill(q.mass());
	    _h_ctheta[2]->fill(ctheta);
	  }
	  else {
	    _h_mpipi[3]->fill(q.mass());
	    _h_ctheta[3]->fill(ctheta);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix) {
	normalize(_h_mpipi [ix]);
	normalize(_h_ctheta[ix]);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_mpipi [4];
    Histo1DPtr _h_ctheta[4];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2017_I1610301);

}
