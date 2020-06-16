// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau -> (5pi)-pi0
  class CLEOII_1994_I373188 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1994_I373188);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_hist,  1, 1, 1);
    }

    void findDecayProducts(const Particle & mother,
                           unsigned int & nstable,
                           Particles& pim, Particles& pi0) {
      for(const Particle & p : mother.children()) {
	long id = p.abspid();
        if (id == PID::PI0 ) {
          pi0.push_back(p);
          ++nstable;
	}
        else if (abs(id) == PID::PIPLUS) {
          pim.push_back(p);
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, pim,pi0);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==PID::TAU)) {
        Particles pi0,pim;
        unsigned int nstable = 0;
        // find the decay products we want
        findDecayProducts(p, nstable, pim, pi0);
        if (nstable != 7) continue;
	// K eta
        if (pim.size() == 5 && pi0.size() == 1) {
	  FourMomentum phad=pi0[0].momentum();
	  for(const Particle & p2: pim)
	    phad += p2.momentum();
	  _hist->fill(phad.mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_hist);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _hist;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1994_I373188);

}
