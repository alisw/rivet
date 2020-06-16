// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau -> pi-pi0 nu_tau
  class CLEO_1999_I508944 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEO_1999_I508944);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_hist_pipi,  1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the taus
      Particles taus;
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if (p.abspid() != PID::TAU) continue;
        Particles pip, pim, pi0;
        unsigned int nstable = 0;
        // get the boost to the rest frame
        // find the decay products we want
        findDecayProducts(p, nstable, pip, pim, pi0);
        if (p.pid() < 0) {
          swap(pip, pim);
        }
        if (nstable != 3) continue;
        // pipi
        if (pim.size() == 1 && pi0.size() == 1)
          _hist_pipi->fill((pi0[0].momentum()+pim[0].momentum()).mass());
      }
    }

    void findDecayProducts(const Particle & mother,
                           unsigned int & nstable,
                           Particles& pip, Particles& pim,
                           Particles& pi0) {
      for(const Particle & p : mother.children()) {
	long id = p.pid();
        if (id == PID::PI0 ) {
          pi0.push_back(p);
          ++nstable;
	}
        else if (id == PID::PIPLUS) {
          pip.push_back(p);
          ++nstable;
        }
        else if (id == PID::PIMINUS) {
          pim.push_back(p);
          ++nstable;
        }
        else if (id == PID::K0S || id == PID::KPLUS ||
		 id == PID::KMINUS) {
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, pip, pim, pi0);
        }
        else
          ++nstable;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_hist_pipi);
      
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist_pipi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CLEO_1999_I508944);


}
