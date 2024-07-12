// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2018_I1679886 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2018_I1679886);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(UnstableParticles(), "UFS");
      book(_h_KK, 1, 1, 1);

    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           unsigned int & nK0, unsigned int & nKp,
			   unsigned int & nKm, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS ) {
	  ++nKp;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::KMINUS ) {
	  ++nKm;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::K0S) {
          ++nK0;
          ++nstable;
	  ptot += p.momentum();
        }
        else if (id == PID::PI0 || id == PID::PIPLUS || id == PID::PIMINUS) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, nK0, nKp, nKm, ptot);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for(const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        unsigned int nstable(0),nK0(0),nKp(0),nKm(0);
      	FourMomentum p_tot(0,0,0,0);
        findDecayProducts(tau, nstable, nK0, nKp, nKm, p_tot);
        if (tau.pid() < 0) {
      	  swap(nKp,nKm);
      	}
       	if(nstable!=3) continue;
      	if(nKm==1 && nK0==1 )
          _h_KK->fill(p_tot.mass());
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_KK;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2018_I1679886);


}
