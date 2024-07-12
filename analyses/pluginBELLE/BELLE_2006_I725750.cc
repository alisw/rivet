// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau -> K phi
  class BELLE_2006_I725750 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I725750);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_hist,  1, 1, 1);
    }

    void findDecayProducts(const Particle & mother,
                           unsigned int & nstable,
                           Particles& phi, Particles& K) {
      for(const Particle & p : mother.children()) {
	long id = p.abspid();
        if (id == PID::PHI ) {
          phi.push_back(p);
          ++nstable;
	}
        else if (id == PID::KPLUS) {
          K.push_back(p);
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, phi,K);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==PID::TAU)) {
        Particles phi, K;
        unsigned int nstable = 0;
        // find the decay products we want
        findDecayProducts(p, nstable, phi, K);
        if (nstable != 3) continue;
	// K phi
        if (K.size() == 1 && phi.size() == 1)
	  _hist->fill((phi[0].momentum()+K[0].momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BELLE_2006_I725750);

}
