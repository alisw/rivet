// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2018_I1697377 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1697377);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_m, 1, 1, 1);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & neta, 
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
        else if (id == PID::ETA) {
	  ++neta;
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, neta,nep,nem,ptot);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over J/psi mesons
      for (const Particle& p :  apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==443)) {
	unsigned nstable(0),neta(0),nep(0),nem(0);
	FourMomentum ptot;
	findDecayProducts(p,nstable,neta,nep,nem,ptot);
	if(nstable==3 && nem==1 && nem==1 && neta==1) {
	  _h_m->fill(ptot.mass());
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_m); // normalize to unity

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_m;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2018_I1697377);


}
