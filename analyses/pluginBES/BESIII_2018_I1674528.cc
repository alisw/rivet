// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Semileptonic Lambda_c+
  class BESIII_2018_I1674528 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1674528);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==4122),"UFS");
      // histos
      book(_h,1,1,1);
    }

    void findDecayProducts(Particle parent, Particles & em, Particles & ep,
			   Particles & nue, Particles & nueBar) {
      for(const Particle & p : parent.children()) {
	if(p.pid() == PID::EMINUS) {
	  em.push_back(p);
	}
	else if(p.pid() == PID::EPLUS) {
	  ep.push_back(p);
	}
	else if(p.pid() == PID::NU_E) {
	  nue.push_back(p);
	}
	else if(p.pid() == PID::NU_EBAR) {
	  nueBar.push_back(p);
	}
	else if(!PID::isHadron(p.pid())) {
	  findDecayProducts(p,em,ep,nue,nueBar);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles()) {
	Particles em,ep,nue,nueBar;
	findDecayProducts(p,em,ep,nue,nueBar);
	if(em.size()==1 && nueBar.size()==1) {
	  double pmod = em[0].momentum().p3().mod();
	  _h->fill(pmod);
	}
	else if(ep.size()==1 && nue.size()==1) {
	  double pmod = ep[0].momentum().p3().mod();
	  _h->fill(pmod);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1674528);

}
