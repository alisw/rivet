// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> KS0 pi+pi+pi-
  class BESIII_2019_I1714778 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1714778);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // histos
      for(unsigned int ix=0;ix<9;++ix)
	book(_h[ix],1,1,1+ix);
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & K0) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PIPLUS) {
	  pip.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KPLUS ||id == PID::KMINUS  ||
		 id == PID::PI0||id == PID::K0L) {
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, K0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid== 411)) {
	unsigned int nstable(0);
	Particles pip, pim, K0;
	findDecayProducts(meson, nstable, pip, pim, K0);
	if (nstable!=4) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if(pip.size()==2&&pim.size()==1&&K0.size()==1) {
	  double mpippim[2] = {(pim[0].momentum()+pip[0].momentum()).mass(), 
			       (pim[0].momentum()+pip[1].momentum()).mass()};
	  if(mpippim[0]>mpippim[1]) {
	    swap(  pip[0],  pip[1]);
	    swap(mpippim[0],mpippim[1]);
	  }
	  _h[0]->fill((K0 [0].momentum()+pim[0].momentum()).mass());
	  _h[1]->fill((K0 [0].momentum()+pip[0].momentum()).mass());
	  _h[2]->fill((K0 [0].momentum()+pip[1].momentum()).mass());
	  _h[3]->fill(mpippim[0]);
	  _h[4]->fill(mpippim[1]);
	  _h[5]->fill((K0 [0].momentum()+pim[0].momentum()+pip[0].momentum()).mass());
	  _h[6]->fill((K0 [0].momentum()+pim[0].momentum()+pip[1].momentum()).mass());
	  _h[7]->fill((pim[0].momentum()+pip[0].momentum()+pip[1].momentum()).mass());
	  _h[8]->fill((pip[0].momentum()+pip[1].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<9;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1714778);

}
