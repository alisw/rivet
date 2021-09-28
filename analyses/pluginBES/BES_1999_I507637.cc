// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> J/Psi pi+pi-
  class BES_1999_I507637 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BES_1999_I507637);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      book(_h_mpipi,1,1,1);
      book(_h_cosl ,1,1,2);
      book(_h_cosX ,1,1,3);
      book(_h_cospi,1,1,4);
    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles & onium) {
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
	else if (id==443) {
	  onium.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,onium);
	}
	else
	  ++nstable;
      }
    }

    void findLeptons(const Particle & mother,
		     unsigned int & nstable,
		     Particles& lp, Particles& lm) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == 11 || id == 13 ) {
	  lm.push_back(p);
	  ++nstable;
	}
       	else if (id == -11 || id==-13) {
       	  lp.push_back(p);
       	  ++nstable;
       	}
	else if ( !p.children().empty() ) {
	  findLeptons(p,nstable,lp,lm);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& psi : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==100443)) {
	unsigned int nstable(0);
	Particles pip, pim, onium;
	findDecayProducts(psi,nstable,pip,pim,onium);
	// check for onium
	if(onium.size() !=1 || nstable !=3 ||
	   pip.size()!=1 || pim.size() !=1 ) continue;
	FourMomentum q = pip[0].momentum()+pim[0].momentum();
	_h_mpipi->fill(q.mass());
	_h_cosX->fill(cos(q.polarAngle()));
	// leptons from J/psi decay
	nstable = 0;
	Particles lp, lm;
	findLeptons(onium[0],nstable,lp,lm);
	if(nstable==2&&lp.size()==1&&lm.size()==1) {
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(onium[0].momentum().betaVec());
	  FourMomentum pl = boost.transform(lp[0].momentum());
	  _h_cosl->fill(cos(pl.polarAngle()));
	}
	// pions in rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(q.betaVec());
	FourMomentum ppi = boost.transform(pip[0].momentum());
	Vector3 axis1 = q.p3().unit();
	_h_cospi->fill(axis1.dot(ppi.p3().unit()));
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mpipi);
      normalize(_h_cosl ,1.,false);
      normalize(_h_cosX );
      normalize(_h_cospi);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_mpipi,_h_cosl,_h_cosX,_h_cospi;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BES_1999_I507637);

}
