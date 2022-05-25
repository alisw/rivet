// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_F1_Decay : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_F1_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");

      // eta pi0 pi0 mode
      book(_h_eta0_etapi0, "eta0_etapi0"   , 70, 0.66, 1.36);
      book(_h_eta0_pi0pi0, "eta0_pi0pi0"   , 80, 0.2, 1.0);
      book(_h_eta0_etapi0pi0, "eta0_etapi0pi0", 70, 1.0, 1.7);
      // eta pi+pi- mode
      book(_h_eta1_etapip, "eta1_etapip"   , 70, 0.66, 1.36);
      book(_h_eta1_etapim, "eta1_etapim"   , 70, 0.66, 1.36);
      book(_h_eta1_pippim, "eta1_pippim"   , 80, 0.2, 1.0);
      book(_h_eta1_etapippim, "eta1_etapippim", 70, 1.0, 1.7);
      // pi+pi-2pi0
      book(_h_4pi0_pi0pi0, "4pi0_pi0pi0"   , 80, 0.2, 1.0);
      book(_h_4pi0_pippi0, "4pi0_pippi0"   , 80, 0.2, 1.0);
      book(_h_4pi0_pimpi0, "4pi0_pimpi0"   , 80, 0.2, 1.0);
      book(_h_4pi0_pippim, "4pi0_pippim"   , 80, 0.2, 1.0);
      book(_h_4pi0_pippimpi0, "4pi0_pippimpi0",100, 0.4, 1.4);
      book(_h_4pi0_pippi0pi0, "4pi0_pippi0pi0",100, 0.4, 1.4);
      book(_h_4pi0_pimpi0pi0, "4pi0_pimpi0pi0",100, 0.4, 1.4);
      book(_h_4pi0_4pi, "4pi0_4pi"      , 70, 1.0, 1.7);
      // 2pi+ 2pi- mode 
      book(_h_4pi1_pippip, "4pi1_pippip"   , 80, 0.2, 1.0);
      book(_h_4pi1_pimpim, "4pi1_pimpim"   , 80, 0.2, 1.0);
      book(_h_4pi1_pippim, "4pi1_pippim"   , 80, 0.2, 1.0);
      book(_h_4pi1_pimpimpip, "4pi1_pimpimpip",100, 0.4, 1.4);
      book(_h_4pi1_pippippim, "4pi1_pippippim",100, 0.4, 1.4);
      book(_h_4pi1_4pi, "4pi1_4pi"      , 70, 1.0, 1.7);
    }

    void findDecayProducts(const Particle & mother,
                           unsigned int & nstable,
                           Particles& pip, Particles& pim,
			   Particles& pi0, Particles& eta) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::ETA ) {
	  eta.push_back(p);
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
	 else if (id == PID::PI0 ) {
           pi0.push_back(p);
	   ++nstable;
	 }
	 else if (id == PID::K0S   || id == PID::K0L ||
		  id == PID::KPLUS || id == PID::KMINUS)
	   ++nstable;
	 else if ( !p.children().empty() ) {
	   findDecayProducts(p,nstable,pip,pim,pi0,eta);
	 }
	 else
	   ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over f_1 mesons
      for(const Particle& f1 : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==20223)) {
	unsigned int nstable(0);
	Particles  pip, pim, pi0, eta;
	findDecayProducts(f1,nstable,pip, pim, pi0, eta);
	// pi+ pi- pi0 pi0
	if(nstable==4 && pip.size()==1 && pim.size()==1 && pi0.size()==2) {
	  _h_4pi0_pi0pi0->fill((pi0[0].momentum()+pi0[1].momentum()).mass(),1.);
	  _h_4pi0_pippi0->fill((pip[0].momentum()+pi0[0].momentum()).mass(),1.);
	  _h_4pi0_pippi0->fill((pip[0].momentum()+pi0[1].momentum()).mass(),1.);
	  _h_4pi0_pimpi0->fill((pim[0].momentum()+pi0[0].momentum()).mass(),1.);
	  _h_4pi0_pimpi0->fill((pim[0].momentum()+pi0[1].momentum()).mass(),1.);
	  _h_4pi0_pippim->fill((pip[0].momentum()+pim[0].momentum()).mass(),1.);	  
	  _h_4pi0_pippimpi0->fill((pip[0].momentum()+pim[0].momentum()+pi0[0].momentum()).mass(),1.);
	  _h_4pi0_pippimpi0->fill((pip[0].momentum()+pim[0].momentum()+pi0[1].momentum()).mass(),1.);
	  _h_4pi0_pippi0pi0->fill((pi0[0].momentum()+pi0[1].momentum()+pip[0].momentum()).mass(),1.);
	  _h_4pi0_pimpi0pi0->fill((pi0[0].momentum()+pi0[1].momentum()+pim[0].momentum()).mass(),1.);
	  _h_4pi0_4pi->fill((pi0[0].momentum()+pi0[1].momentum()+pim[0].momentum()+pip[0].momentum()).mass(),1.);
	}
	else if(nstable==4 && pip.size()==2 && pim.size()==2) {
	  _h_4pi1_pippip   ->fill((pip[0].momentum()+pip[1].momentum()).mass(),1.);
	  _h_4pi1_pimpim   ->fill((pim[0].momentum()+pim[1].momentum()).mass(),1.);
	  _h_4pi1_pippim   ->fill((pip[0].momentum()+pim[0].momentum()).mass(),1.);
	  _h_4pi1_pippim   ->fill((pip[0].momentum()+pim[1].momentum()).mass(),1.);
	  _h_4pi1_pippim   ->fill((pip[1].momentum()+pim[0].momentum()).mass(),1.);
	  _h_4pi1_pippim   ->fill((pip[1].momentum()+pim[1].momentum()).mass(),1.);
	  _h_4pi1_pimpimpip->fill((pim[0].momentum()+pim[1].momentum()+pip[0].momentum()).mass(),1.);
	  _h_4pi1_pimpimpip->fill((pim[0].momentum()+pim[1].momentum()+pip[1].momentum()).mass(),1.);
	  _h_4pi1_pippippim->fill((pip[0].momentum()+pip[1].momentum()+pim[0].momentum()).mass(),1.);
	  _h_4pi1_pippippim->fill((pip[0].momentum()+pip[1].momentum()+pim[1].momentum()).mass(),1.);
	  _h_4pi1_4pi      ->fill((pip[0].momentum()+pip[1].momentum()+
				   pim[0].momentum()+pim[1].momentum()).mass(),1.);
	}
	else if(nstable==3 && eta.size()==1 && pip.size()==1 && pim.size()==1) {
	  _h_eta1_etapip   ->fill((eta[0].momentum()+pip[0].momentum()).mass(),1.);
	  _h_eta1_etapim   ->fill((eta[0].momentum()+pim[0].momentum()).mass(),1.);
	  _h_eta1_pippim   ->fill((pim[0].momentum()+pip[0].momentum()).mass(),1.);
	  _h_eta1_etapippim->fill((eta[0].momentum()+pim[0].momentum()+pip[0].momentum()).mass(),1.);
	}
	else if(nstable==3 && eta.size()==1 && pi0.size()==2 ) {
	  _h_eta0_etapi0   ->fill((eta[0].momentum()+pi0[0].momentum()).mass(),1.);
	  _h_eta0_etapi0   ->fill((eta[0].momentum()+pi0[1].momentum()).mass(),1.);
	  _h_eta0_pi0pi0   ->fill((pi0[0].momentum()+pi0[1].momentum()).mass(),1.);
	  _h_eta0_etapi0pi0->fill((eta[0].momentum()+pi0[0].momentum()+pi0[1].momentum()).mass(),1.);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_eta0_etapi0   );
      normalize(_h_eta0_pi0pi0   );
      normalize(_h_eta0_etapi0pi0);

      normalize(_h_eta1_etapip);
      normalize(_h_eta1_etapim);
      normalize(_h_eta1_pippim);
      normalize(_h_eta1_etapippim);

      normalize(_h_4pi0_pi0pi0);
      normalize(_h_4pi0_pippi0);
      normalize(_h_4pi0_pimpi0);
      normalize(_h_4pi0_pippim);
      normalize(_h_4pi0_pippimpi0);
      normalize(_h_4pi0_pippi0pi0);
      normalize(_h_4pi0_pimpi0pi0);
      normalize(_h_4pi0_4pi);
      
      normalize(_h_4pi1_pippip);
      normalize(_h_4pi1_pimpim);
      normalize(_h_4pi1_pippim);
      normalize(_h_4pi1_pimpimpip);
      normalize(_h_4pi1_pippippim);
      normalize(_h_4pi1_4pi);
    }

    //@}


    // @name Histograms
    //@{
    Histo1DPtr _h_eta0_etapi0,_h_eta0_pi0pi0,_h_eta0_etapi0pi0;
    Histo1DPtr _h_eta1_etapip,_h_eta1_etapim,_h_eta1_pippim,_h_eta1_etapippim;
    Histo1DPtr _h_4pi0_pi0pi0,_h_4pi0_pippi0,_h_4pi0_pimpi0,_h_4pi0_pippim,_h_4pi0_pippimpi0,_h_4pi0_pippi0pi0,
      _h_4pi0_pimpi0pi0,_h_4pi0_4pi;
    Histo1DPtr _h_4pi1_pippip,_h_4pi1_pimpim,_h_4pi1_pippim,_h_4pi1_pimpimpip,_h_4pi1_pippippim,_h_4pi1_4pi;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_F1_Decay);


}
