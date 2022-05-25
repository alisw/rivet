// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2007_I753243 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I753243);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(UnstableParticles(), "UFS");
      book(_hist_Kpi, 1, 1, 1);

    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           unsigned int & npip, unsigned int & npim,
                           unsigned int & nK, FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::K0S ) {
	  ++nK;
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
        else if (id == PID::PI0 || id == PID::KPLUS ||
		 id == PID::KMINUS) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, npip, npim, nK,ptot);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for(const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        unsigned int nstable(0),npip(0),npim(0),nK(0);
	FourMomentum p_tot(0,0,0,0);
        findDecayProducts(tau, nstable, npip, npim, nK, p_tot);
        if (tau.pid() < 0) swap(npip, npim);
	if(nstable==3 && npim==1 && nK==1)
          _hist_Kpi->fill(p_tot.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_hist_Kpi);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist_Kpi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BELLE_2007_I753243);


}
