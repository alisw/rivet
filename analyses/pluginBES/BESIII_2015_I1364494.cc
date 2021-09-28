// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2015_I1364494 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2015_I1364494);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_m, 1, 1, 3);
      book(_netap, "TMP/netap");
    }
    
    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & ngamma, 
                           unsigned int & nep, unsigned int & nem, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::EMINUS ) {
	  ++nem;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::EPLUS) {
          ++nep;
          ++nstable;
	  ptot += p.momentum();
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, ngamma,nep,nem,ptot);
        }
        else if (id == PID::GAMMA) {
	  ++ngamma;
          ++nstable;
        }
        else
          ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over eta' mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==331)) {
	unsigned nstable(0),ngamma(0),nep(0),nem(0);
	FourMomentum ptot;
	findDecayProducts(p,nstable,ngamma,nep,nem,ptot);
	if(nstable==3 && nem==1 && nem==1 && ngamma==1)
	  _h_m->fill(ptot.mass());
	else if(nstable==2 &&ngamma==2)
	  _netap->fill();
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // divide by no so BR and mult by bin width
      // and 100 as in %
      scale(_h_m,1./_netap->sumW()*0.1*100.);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_m;
    CounterPtr _netap;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESIII_2015_I1364494);


}
