// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D0 -> eta pi+pi- decays
  class CLEOC_2008_I779705 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2008_I779705);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(),"UFS");
      // histograms
      book(_h_eta_pi,1,1,1);
      book(_h_pi_pi ,1,1,2);
    }

    void findDecay(const Particle & parent, unsigned int & nstable,
		   Particles & pip, Particles & pim, Particles & eta) {
      for(const Particle & p : parent.children()) {
	if(p.pid()==PID::PIPLUS) {
	  ++nstable;
	  pip.push_back(p);
	}
	else if(p.pid()==PID::PIMINUS) {
	  ++nstable;
	  pim.push_back(p);
	}
	else if(p.pid()==PID::ETA) {
	  ++nstable;
	  eta.push_back(p);
	}
	else if(p.children().empty() || p.pid()==PID::K0S) {
	  ++nstable;
	}
	else {
	  findDecay(p,nstable,pip,pim,eta);
	}	
      }
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & D0 : apply<UnstableParticles>(event,"UFS").particles(Cuts::abspid==421)) {
	unsigned int nstable(0);
	Particles pip, pim, eta;
	findDecay(D0,nstable,pip,pim,eta);
	if(nstable==3 && pip.size()==1 && pim.size()==1 && eta.size()==1) {
	  if(D0.pid()<0) swap (pip,pim);
	  _h_eta_pi->fill((pip[0].momentum()+eta[0].momentum()).mass());
	  _h_pi_pi ->fill((pip[0].momentum()+pim[0].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_eta_pi);
      normalize(_h_pi_pi );
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_eta_pi,_h_pi_pi;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2008_I779705);

}
