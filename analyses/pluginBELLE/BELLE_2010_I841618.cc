// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2010_I841618 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2010_I841618);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(UnstableParticles(), "UFS");
      book(_h_3pi,   1, 1, 1);
      book(_h_Kpipi, 2, 1, 1);
      book(_h_KKpi,  3, 1, 1);
      book(_h_3K,    4, 1, 1);

    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           unsigned int & npip, unsigned int & npim,
                           unsigned int & nKp, unsigned int & nKm, FourMomentum & ptot) {
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
        else if (id == PID::PIPLUS) {
          ++npip;
          ++nstable;
	  ptot += p.momentum();
        }
        else if (id == PID::PIMINUS) {
          ++npim;
          ++nstable;
	  ptot += p.momentum();
        }
        else if (id == PID::PI0 || id == PID::K0S) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, npip, npim, nKp, nKm,ptot);
        }
        else
          ++nstable;
      }
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for(const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        unsigned int nstable(0),npip(0),npim(0),nKp(0),nKm(0);
      	FourMomentum p_tot(0,0,0,0);
        findDecayProducts(tau, nstable, npip, npim, nKp, nKm, p_tot);
        if (tau.pid() < 0) {
      	  swap(npip, npim);
      	  swap(nKp,nKm);
      	}
       	if(nstable!=4) continue;
      	if(npim==2 && npip==1 )
          _h_3pi->fill(p_tot.mass());
      	else if(npim==1 && npip==1 && nKm==1)
          _h_Kpipi->fill(p_tot.mass());
      	else if(nKm==1 && nKp==1 && npim==1)
          _h_KKpi->fill(p_tot.mass());
      	else if(nKm==2 && nKp==1 )
          _h_3K->fill(p_tot.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_3pi); 
      normalize(_h_Kpipi); 
      normalize(_h_KKpi); 
      normalize(_h_3K); 

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_3pi, _h_Kpipi, _h_KKpi, _h_3K;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BELLE_2010_I841618);


}
