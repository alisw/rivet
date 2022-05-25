// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief phi -> eta pi0 gamma
  class KLOE_2009_I818106 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(KLOE_2009_I818106);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_h_etapi,  1, 1, 1);
      book(_nPhi, "TMP/PhiCounter");
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           unsigned int & neta,   unsigned int & npi,
			   unsigned int & ngamma, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::ETA ) {
	  ++neta;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::PI0) {
	  ++npi;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::GAMMA) {
	  ++ngamma;
          ++nstable;
	}
        else if (id == PID::PIPLUS || id == PID::PIMINUS) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, neta, npi, ngamma, ptot);
        }
        else
          ++nstable;
      }
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over phis
      for(const Particle& phi : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::PHI)) {
	_nPhi->fill();
        unsigned int nstable(0),neta(0),npi(0),ngamma(0);
      	FourMomentum p_tot(0,0,0,0);
        findDecayProducts(phi, nstable, neta, npi, ngamma, p_tot);
       	if(nstable!=3) continue;
      	if(neta==1 && npi==1 && ngamma==1 ) {
          _h_etapi->fill(p_tot.mass()/MeV);
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalise to total no of phi mesons
      // and mult by 10^7 due normalisation in paper
      scale( _h_etapi, 1./_nPhi->sumW()*1e7);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_etapi;
    CounterPtr _nPhi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(KLOE_2009_I818106);


}
