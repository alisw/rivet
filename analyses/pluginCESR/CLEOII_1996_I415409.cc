// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau -> K eta
  class CLEOII_1996_I415409 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1996_I415409);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_hist,  1, 1, 1);
    }

    void findDecayProducts(const Particle & mother,
                           unsigned int & nstable,
                           Particles& eta, Particles& K) {
      for(const Particle & p : mother.children()) {
	long id = p.abspid();
        if (id == PID::ETA ) {
          eta.push_back(p);
          ++nstable;
	}
        else if (id == PID::KPLUS) {
          K.push_back(p);
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, eta,K);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==PID::TAU)) {
        Particles eta, K;
        unsigned int nstable = 0;
        // find the decay products we want
        findDecayProducts(p, nstable, eta, K);
        if (nstable != 3) continue;
	// K eta
        if (K.size() == 1 && eta.size() == 1)
	  _hist->fill((eta[0].momentum()+K[0].momentum()).mass());
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


  DECLARE_RIVET_PLUGIN(CLEOII_1996_I415409);

}
