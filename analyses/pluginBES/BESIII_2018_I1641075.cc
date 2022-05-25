// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2018_I1641075 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1641075);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_m, 1, 1, 5);

    }
    
    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & ngamma, 
                           unsigned int & npip, unsigned int & npim, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PIMINUS ) {
	  ++npim;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::PIPLUS) {
          ++npip;
          ++nstable;
	  ptot += p.momentum();
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, ngamma,npip,npim,ptot);
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
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==331)) {
	unsigned nstable(0),ngamma(0),npip(0),npim(0);
	FourMomentum ptot;
	findDecayProducts(p,nstable,ngamma,npip,npim,ptot);
	if(nstable==3 && npim==1 && npip==1 && ngamma==1)
	  _h_m->fill(ptot.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_m);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_m;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2018_I1641075);


}
