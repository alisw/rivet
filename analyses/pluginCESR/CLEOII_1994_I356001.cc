// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CLEOII_1994_I356001 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1994_I356001);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      // Histograms
      book(_n_2S,"/TMP/n_2S");
      book(_n_3S,"/TMP/n_3S");
      book(_h_3S_1S_neutral,1,1,1);
      book(_h_3S_2S_neutral,1,1,2);
      book(_h_3S_1S_ex     ,2,1,1);
      book(_h_3S_1S_in     ,2,1,2);
      book(_h_3S_2S_ex     ,3,1,1);
      book(_h_3S_2S_in     ,3,1,2);
      book(_h_2S_1S_ex     ,4,1,1);
      book(_h_2S_1S_in     ,4,1,2);
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
	else if (abs(id)%1000==553) {
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
      // loop over unstable particles
      for(const Particle& ups : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==100553 or
										 Cuts::pid==200553)) {
	if(ups.pid()==200553)
	  _n_3S->fill();
	else
	  _n_2S->fill();
	unsigned int nstable(0);
	Particles pip, pim, pi0, onium;
	findDecayProducts(ups,nstable,pip,pim,pi0,onium);
	// check for onium
	if(onium.size() !=1 || nstable !=3) continue;
	//pi+pi-
	if( pip.size()==1 && pim.size() ==1) {
	  double q = (pip[0].momentum()+pim[0].momentum()).mass();
	  if(ups.pid()==200553 && onium[0].pid()==553) {
	    _h_3S_1S_ex->fill(q);
	    _h_3S_1S_in->fill(q);
	  }
	  else if(ups.pid()==200553 && onium[0].pid()==100553) {
	    _h_3S_2S_ex->fill(q);
	    _h_3S_2S_in->fill(q);
	  }
	  else if(ups.pid()==100553 && onium[0].pid()==553) {
	    _h_2S_1S_ex->fill(q);
	    _h_2S_1S_in->fill(q);
	  }
	}
	else if(pi0.size()==2) {
	  double q = (pi0[0].momentum()+pi0[1].momentum()).mass();
	  if(ups.pid()==200553 && onium[0].pid()==553) {
	    _h_3S_1S_neutral->fill(q);
	  }
	  else if(ups.pid()==200553 && onium[0].pid()==100553) {
	    _h_3S_2S_neutral->fill(q);
	  }

	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // widths of Upsilon 2S and 3S
      vector<double> width = {31.98,20.32};
      scale(_h_3S_2S_neutral, width[1]/ *_n_3S);
      scale(_h_3S_1S_neutral, width[1]/ *_n_3S);
      scale(_h_3S_2S_ex     , width[1]/ *_n_3S);
      scale(_h_3S_2S_in     , width[1]/ *_n_3S);
      scale(_h_3S_1S_ex     , width[1]/ *_n_3S);
      scale(_h_3S_1S_in     , width[1]/ *_n_3S);
      scale(_h_2S_1S_ex     , width[0]/ *_n_2S);
      scale(_h_2S_1S_in     , width[0]/ *_n_2S);
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _n_3S,_n_2S;
    Histo1DPtr _h_3S_2S_neutral,_h_3S_1S_neutral;
    Histo1DPtr _h_3S_2S_ex,_h_3S_2S_in,_h_3S_1S_ex,_h_3S_1S_in;
    Histo1DPtr _h_2S_1S_ex,_h_2S_1S_in;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1994_I356001);

}
