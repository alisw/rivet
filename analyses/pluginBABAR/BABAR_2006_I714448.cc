// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2006_I714448 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2006_I714448);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      book(_h_1S, 1, 1, 1);
      book(_h_2S, 1, 1, 2);
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
	else if (abs(id)%1000==553) {
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

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==300553)) {
	unsigned int nstable(0);
	Particles pip, pim, onium;
	findDecayProducts(ups,nstable,pip,pim,onium);
	// check for onium
	if(onium.size() !=1 || nstable !=3) continue;
	// check for pipi
	if( ! (pip.size()==1 && pim.size() ==1) ) continue;
	FourMomentum q = pip[0].momentum()+pim[0].momentum();
	if(onium[0].pid()==553)
	  _h_1S->fill(q.mass());
	else if(onium[0].pid()==100553)
	  _h_2S->fill(q.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_1S);
      normalize(_h_2S);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_1S,_h_2S;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2006_I714448);

}
