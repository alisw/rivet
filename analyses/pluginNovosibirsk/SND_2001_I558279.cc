// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief phi -> 3 pion decay
  class SND_2001_I558279 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2001_I558279);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_pm, 1, 1, 1);
      book(_h_p0, 2, 1, 1);
      
    }
    
    void findDecayProducts(const Particle & mother, unsigned int & nstable, Particles &pip,
			   Particles &pim, Particles &pi0) {
      for(const Particle & p : mother.children()) {
	int id = p.pid();
	if (id == PID::PIMINUS ) {
	  pim.push_back(p);
	  ++nstable;
	}
        else if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, pip,pim,pi0);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==333)) {
	Particles pip,pim,pi0;
	unsigned int nstable(0);
	findDecayProducts(p, nstable, pip,pim,pi0);
	if(nstable==3 && pip.size()==1 && pim.size()==1&&pi0.size()==1) {
	  _h_pm->fill((pip[0].momentum()+pim[0].momentum()).mass()/MeV);
	  _h_p0->fill((pip[0].momentum()+pi0[0].momentum()).mass()/MeV);
	  _h_p0->fill((pim[0].momentum()+pi0[0].momentum()).mass()/MeV);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_pm);
      normalize(_h_p0);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_pm,_h_p0;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(SND_2001_I558279);


}
