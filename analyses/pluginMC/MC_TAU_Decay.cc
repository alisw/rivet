// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_TAU_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_TAU_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(),"UFS");

      // Book histograms
      // leptonic
      book(_h_2B_m2enu, "h_2B_m2enu", 200,0.,3.15);
      book(_h_2B_menu, "h_2B_menu" , 200,0.,1.8 );
      // hadronic
      // 1 hadron
      book(_h_1B_xpi, "h_1B_xpi", 25,0.0,1.0);
      // 2 hadrons
      book(_h_2B_m2pipi, "h_2B_m2pipi", 200,0.,3.15);
      book(_h_2B_mpipi, "h_2B_mpipi" , 200,0.,1.8 );
      book(_h_2B_m2munu, "h_2B_m2munu", 200,0.,3.15);
      book(_h_2B_mmunu, "h_2B_mmunu" , 200,0.,1.8 );
      book(_h_2B_m2KpiA, "h_2B_m2KpiA", 200,0.,3.15);
      book(_h_2B_mKpiA, "h_2B_mKpiA" , 200,0.,1.8 );
      book(_h_2B_m2KpiB, "h_2B_m2KpiB", 200,0.,3.15);
      book(_h_2B_mKpiB, "h_2B_mKpiB" , 200,0.,1.8 );
      book(_h_2B_m2Keta, "h_2B_m2Keta", 200,0.,3.15);
      book(_h_2B_mKeta, "h_2B_mKeta" , 200,0.,1.8 );
      book(_h_2B_m2KK, "h_2B_m2KK"  , 200,0.,3.15);
      book(_h_2B_mKK, "h_2B_mKK"   , 200,0.,1.8 );
      // 3 hadrons
      Histo1DPtr dummy;
      for(unsigned int ix=0;ix<4;++ix) {
	if(ix<3) {
	  std::ostringstream title1 ; title1  << "h_3B_pippimpim_" << ix+1;
          book(dummy, title1 .str(),200,0.,1.8);
	  _h_3B_pippimpim   .push_back(dummy);
	  std::ostringstream title2 ; title2  << "h_3B_pi0pi0pim_" << ix+1;
          book(dummy, title2 .str(),200,0.,1.8);
	  _h_3B_pi0pi0pim   .push_back(dummy);
	  std::ostringstream title5 ; title5  << "h_3B_pi0pi0km_" << ix+1;
          book(dummy, title5 .str(),200,0.,1.8);
	  _h_3B_pi0pi0km    .push_back(dummy);
	  std::ostringstream title10; title10 << "h_3B_kspimks_" << ix+1;
          book(dummy, title10.str(),200,0.,1.8);
	  _h_3B_kspimks     .push_back(dummy);
	  std::ostringstream title11; title11 << "h_3B_klpimkl_" << ix+1;
          book(dummy, title11.str(),200,0.,1.8);
	  _h_3B_klpimkl     .push_back(dummy);
	}
	std::ostringstream title3 ; title3  << "h_3B_kmpimkp_" << ix+1;
        book(dummy, title3.str(),200,0.,1.8);
	_h_3B_kmpimkp     .push_back(dummy);
	std::ostringstream title4 ; title4  << "h_3B_kmpi0k0_" << ix+1;
        book(dummy, title4.str(),200,0.,1.8);
	_h_3B_kmpi0k0     .push_back(dummy);
	std::ostringstream title6 ; title6  << "h_3B_kmpimpip_" << ix+1;
        book(dummy, title6.str(),200,0.,1.8);
	_h_3B_kmpimpip    .push_back(dummy);
	std::ostringstream title7 ; title7  << "h_3B_pimk0pi0_" << ix+1;
        book(dummy, title7.str(),200,0.,1.8);
	_h_3B_pimk0pi0    .push_back(dummy);
	std::ostringstream title8 ; title8  << "h_3B_pimpi0eta_" << ix+1;
        book(dummy, title8.str(),200,0.,1.8);
	_h_3B_pimpi0eta   .push_back(dummy);
	std::ostringstream title9 ; title9  << "h_3B_pimpi0gamma_" << ix+1;
        book(dummy, title9.str(),200,0.,1.8);
	_h_3B_pimpi0gamma .push_back(dummy);
	std::ostringstream title12; title12 << "h_3B_kspimkl_" << ix+1;
        book(dummy, title12.str(),200,0.,1.8);
	_h_3B_kspimkl     .push_back(dummy);
      }
      // 4 pion decays
      for(unsigned int ix=0;ix<5;++ix) {
	std::ostringstream title1 ; title1  << "h_4B_pipi_" << ix+1;
        book(dummy, title1.str(),200,0.,1.8);
	_h_4B_pipi  .push_back(dummy);
	std::ostringstream title2 ; title2  << "h_4B_pipipi_" << ix+1;
        book(dummy, title2.str(),200,0.,1.8);
	_h_4B_pipipi.push_back(dummy);
      }
      book(dummy, "h_4B_pipi_6",200,0.,1.8);
      _h_4B_pipi  .push_back(dummy);
      for(unsigned int ix=0;ix<2;++ix) {
	std::ostringstream title ; title  << "h_4B_pipipipi_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_4B_pipipipi.push_back(dummy);
      }
      // 5 pion decays
      // 2 pi0 2pi- pi+
      book(_h_5B_q1, "h_5B_q1",200,0.,1.8);
      for(unsigned int ix=0;ix<5;++ix) {
	std::ostringstream title;
	title << "h_5B_pipi1_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipi1.push_back(dummy);
      }
      for(unsigned int ix=0;ix<5;++ix) {
	std::ostringstream title;
	title << "h_5B_pipipi1_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipipi1.push_back(dummy);
      }
      for(unsigned int ix=0;ix<3;++ix) {
	std::ostringstream title;
	title << "h_5B_pipipipi1_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipipipi1.push_back(dummy);
      }
      // 4 pi0  pi-
      book(_h_5B_q2, "h_5B_q2",200,0.,1.8);
      for(unsigned int ix=0;ix<2;++ix) {
	std::ostringstream title;
	title << "h_5B_pipi2_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipi2.push_back(dummy);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	std::ostringstream title;
	title << "h_5B_pipipi2_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipipi2.push_back(dummy);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	std::ostringstream title;
	title << "h_5B_pipipipi2_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipipipi2.push_back(dummy);
      }
      // 3 pi- 2 pi+
      book(_h_5B_q3, "h_5B_q3",200,0.,1.8);
      for(unsigned int ix=0;ix<3;++ix) {
	std::ostringstream title;
	title << "h_5B_pipi3_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipi3.push_back(dummy);
      }
      for(unsigned int ix=0;ix<3;++ix) {
	std::ostringstream title;
	title << "h_5B_pipipi3_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipipi3.push_back(dummy);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	std::ostringstream title;
	title << "h_5B_pipipipi3_" << ix+1;
        book(dummy, title.str(),200,0.,1.8);
	_h_5B_pipipipi3.push_back(dummy);
      }
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & ep  , Particles & em  , Particles & nu_e , Particles & nu_ebar,
			   Particles & mup , Particles & mum , Particles & nu_mu, Particles & nu_mubar,
			   Particles & pip , Particles & pim , Particles & pi0  ,
			   Particles & Kp  , Particles & Km  , Particles & K0S  , Particles & K0L,
			   Particles & eta , Particles & gamma) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if ( id == PID::KPLUS ) {
       	  Kp.push_back(p);
	  ++nstable;
	}
	else if (id == PID::KMINUS ) {
	  Km.push_back(p);
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
	else if (id == PID::EPLUS) {
	  ep.push_back(p);
	  ++nstable;
	}
	else if (id == PID::EMINUS) {
	  em.push_back(p);
	  ++nstable;
	}
	else if (id == PID::NU_E) {
	  nu_e.push_back(p);
	  ++nstable;
	}
	else if (id == PID::NU_EBAR) {
	  nu_ebar.push_back(p);
	  ++nstable;
	}
	else if (id == PID::NU_MU) {
	  nu_mu.push_back(p);
	  ++nstable;
	}
	else if (id == PID::NU_MUBAR) {
	  nu_mubar.push_back(p);
	  ++nstable;
	}
	else if (id == PID::ANTIMUON) {
	  mup.push_back(p);
	  ++nstable;
	}
	else if (id == PID::MUON) {
	  mum.push_back(p);
	  ++nstable;
	}
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if (id == PID::K0S) {
	  K0S.push_back(p);
          ++nstable;
        }
	else if (id == PID::K0L) {
	  K0L.push_back(p);
          ++nstable;
        }
	else if (id == PID::ETA) {
	  eta.push_back(p);
          ++nstable;
        }
	else if (id == PID::PHOTON) {
	  gamma.push_back(p);
          ++nstable;
        }
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable,ep,em,nu_e,nu_ebar,mup,mum,nu_mu,nu_mubar,
			    pip, pim, pi0,Kp , Km, K0S, K0L,eta,gamma);
	}
	else
	  ++nstable;
      }
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& tau : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
	unsigned int nstable(0);
	Particles ep,em,nu_e,nu_ebar,mup,mum,nu_mu,nu_mubar;
	Particles pip, pim, pi0, Kp , Km, K0S, K0L, eta,gamma;
	findDecayProducts(tau, nstable,ep,em,nu_e,nu_ebar,mup,mum,nu_mu,nu_mubar,
			  pip, pim, pi0,Kp , Km, K0S, K0L,eta,gamma);
	if(tau.pid()<0) {
	  swap(pim,pip);
	  swap(Kp,Km);
	  swap(em,ep);
	  swap(mum,mup);
	  swap(nu_e ,nu_ebar );
	  swap(nu_mu,nu_mubar);
	}
	// cerr << "testing before loop " << nstable << " "
	//      << pip.size() << " " << pim.size() << " " << pi0.size() << " " 
	//      << Kp.size() << " " << Km.size() << " " << K0S.size() << " " <<  K0L.size() << "\n";
	// 2 hadrons
        if(nstable==2) {
          if(pim.size()==1) {
            double xpi = pim[0].momentum().p()/tau.momentum().p();
            _h_1B_xpi->fill(xpi);
          }
        }
	else if(nstable==3 ) {
	  if(em.size()==1 && nu_ebar.size()==1) {
	    FourMomentum ptot = em[0].momentum()+nu_ebar[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2enu->fill(     mass2 );
	    _h_2B_menu ->fill(sqrt(mass2));
	  }
	  else if(mum.size()==1 && nu_mubar.size()==1) {
	    FourMomentum ptot = mum[0].momentum()+nu_mubar[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2munu->fill(     mass2 );
	    _h_2B_mmunu ->fill(sqrt(mass2));
	  }
	  else if(pim.size()==1 && pi0.size()==1) {
	    FourMomentum ptot = pim[0].momentum()+pi0[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2pipi->fill(     mass2 );
	    _h_2B_mpipi ->fill(sqrt(mass2));
	  }
	  else if(Km.size()==1 && pi0.size()==1) {
	    FourMomentum ptot = Km[0].momentum()+pi0[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2KpiA->fill(     mass2 );
	    _h_2B_mKpiA ->fill(sqrt(mass2));
	  }
	  else if(K0S.size()==1&&pim.size()==1) {
	    FourMomentum ptot = K0S[0].momentum()+pim[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2KpiB->fill(     mass2 );
	    _h_2B_mKpiB ->fill(sqrt(mass2));
	  }
	  else if(K0L.size()==1&&pim.size()==1) {
	    FourMomentum ptot = K0L[0].momentum()+pim[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2KpiB->fill(     mass2 );
	    _h_2B_mKpiB ->fill(sqrt(mass2));
	  }
	  else if(K0S.size()==1&&Km.size()==1) {
	    FourMomentum ptot = K0S[0].momentum()+Km[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2KK->fill(     mass2 );
	    _h_2B_mKK ->fill(sqrt(mass2));
	  }
	  else if(K0L.size()==1&&Km.size()==1) {
	    FourMomentum ptot = K0L[0].momentum()+Km[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2KK->fill(     mass2 );
	    _h_2B_mKK ->fill(sqrt(mass2));
	  }
	  else if(eta.size()==1&&Km.size()==1) {
	    FourMomentum ptot = eta[0].momentum()+Km[0].momentum();
	    double mass2 = ptot.mass2();
	    _h_2B_m2Keta->fill(     mass2 );
	    _h_2B_mKeta ->fill(sqrt(mass2));
	  }
	}
	else if(nstable==4) {
	  if(pim.size()==2&&pip.size()==1) {
	    _h_3B_pippimpim[0]->fill((pim[0].momentum()+pim[1].momentum()+pip[0].momentum()).mass());
	    _h_3B_pippimpim[1]->fill((pim[0].momentum()+pim[1].momentum()).mass());
	    _h_3B_pippimpim[2]->fill((pim[0].momentum()+pip[0].momentum()).mass());
	    _h_3B_pippimpim[2]->fill((pim[1].momentum()+pip[0].momentum()).mass());
	  }
	  else if(pim.size()==1&&pi0.size()==2) {
	    _h_3B_pi0pi0pim[0]->fill((pi0[0].momentum()+pi0[1].momentum()+pim[0].momentum()).mass());
	    _h_3B_pi0pi0pim[1]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	    _h_3B_pi0pi0pim[2]->fill((pi0[0].momentum()+pim[0].momentum()).mass());
	    _h_3B_pi0pi0pim[2]->fill((pi0[1].momentum()+pim[0].momentum()).mass());
	  }
	  else if(Km.size()==1&&Kp.size()==1&&pim.size()==1) {
	    _h_3B_kmpimkp[0]->fill((Km[0].momentum()+pim[0].momentum()+Kp[0].momentum()).mass());
	    _h_3B_kmpimkp[1]->fill((Km[0].momentum()+pim[0].momentum()).mass());
	    _h_3B_kmpimkp[2]->fill((Km[0].momentum()+ Kp[0].momentum()).mass());
	    _h_3B_kmpimkp[3]->fill((Kp[0].momentum()+pim[0].momentum()).mass());
	  }
	  else if((K0S.size()==1||K0L.size()==1)&&Km.size()==1&&pi0.size()==1) {
	    FourMomentum pk = K0L.size()==1 ? K0L[0].momentum() : K0S[0].momentum();
	    _h_3B_kmpi0k0[0]->fill((Km[0].momentum()+pi0[0].momentum()+pk).mass());
	    _h_3B_kmpi0k0[1]->fill((Km[0].momentum()+pi0[0].momentum()).mass());
	    _h_3B_kmpi0k0[2]->fill((Km[0].momentum()+pk ).mass());
	    _h_3B_kmpi0k0[3]->fill((pk+pi0[0].momentum()).mass());
	  }
	  else if(pi0.size()==2&&Km.size()==1) {
	    _h_3B_pi0pi0km[0]->fill((pi0[0].momentum()+pi0[1].momentum()+Km[0].momentum()).mass());
	    _h_3B_pi0pi0km[1]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	    _h_3B_pi0pi0km[2]->fill((pi0[0].momentum()+Km[0].momentum() ).mass());
	    _h_3B_pi0pi0km[2]->fill((pi0[1].momentum()+Km[0].momentum() ).mass());
	  }
	  else if(Km.size()==1&&pim.size()==1&&pip.size()==1) {
	    _h_3B_kmpimpip[0]->fill((pip[0].momentum()+pim[0].momentum()+Km[0].momentum()).mass());
	    _h_3B_kmpimpip[1]->fill((Km[0].momentum()+pim[0].momentum()).mass());
	    _h_3B_kmpimpip[2]->fill((Km[0].momentum()+pip[0].momentum() ).mass());
	    _h_3B_kmpimpip[3]->fill((pip[0].momentum()+pim[0].momentum() ).mass());
	  }
	  else if(pim.size()==1&&(K0S.size()==1||K0L.size()==1)&&pi0.size()==1) {
	    FourMomentum pk = K0L.size()==1 ? K0L[0].momentum() : K0S[0].momentum();
	    _h_3B_pimk0pi0[0]->fill((pim[0].momentum()+pi0[0].momentum()+pk).mass());
	    _h_3B_pimk0pi0[1]->fill((pim[0].momentum()+pk).mass());
	    _h_3B_pimk0pi0[2]->fill((pim[0].momentum()+pi0[0].momentum()  ).mass());
	    _h_3B_pimk0pi0[3]->fill((pk+pi0[0].momentum()).mass());
	  }
	  else if(pim.size()==1&&pi0.size()==1&&eta.size()==1) {
	    _h_3B_pimpi0eta[0]->fill((pim[0].momentum()+pi0[0].momentum()+eta[0].momentum()).mass());
	    _h_3B_pimpi0eta[1]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	    _h_3B_pimpi0eta[2]->fill((pim[0].momentum()+eta[0].momentum()).mass());
	    _h_3B_pimpi0eta[3]->fill((pi0[0].momentum()+eta[0].momentum()).mass());
	  }
	  else if(pim.size()==1&&pi0.size()==1&&gamma.size()==1) {
	    _h_3B_pimpi0gamma[0]->fill((pim[0].momentum()+pi0[0].momentum()+gamma[0].momentum()).mass());
	    _h_3B_pimpi0gamma[1]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	    _h_3B_pimpi0gamma[2]->fill((pim[0].momentum()+gamma[0].momentum()).mass());
	    _h_3B_pimpi0gamma[3]->fill((pi0[0].momentum()+gamma[0].momentum()).mass());
	  }
	  else if(K0S.size()==2&&pim.size()==1) {
	    _h_3B_kspimks[0]->fill((pim[0].momentum()+K0S[0].momentum()+K0S[1].momentum()).mass());
	    _h_3B_kspimks[1]->fill((pim[0].momentum()+K0S[0].momentum()).mass());
	    _h_3B_kspimks[1]->fill((pim[0].momentum()+K0S[1].momentum()).mass());
	    _h_3B_kspimks[2]->fill((K0S [0].momentum()+K0S[1].momentum()).mass());
	  }
	  else if(K0L.size()==2&&pim.size()==1) {
	    _h_3B_klpimkl[0]->fill((pim[0].momentum()+K0L[0].momentum()+K0L[1].momentum()).mass());
	    _h_3B_klpimkl[1]->fill((pim[0].momentum()+K0L[0].momentum()).mass());
	    _h_3B_klpimkl[1]->fill((pim[0].momentum()+K0L[1].momentum()).mass());
	    _h_3B_klpimkl[2]->fill((K0L [0].momentum()+K0L[1].momentum()).mass());
	  }
	  else if(K0S.size()==1&&K0L.size()==1&&pim.size()==1) {
	    _h_3B_kspimkl[0]->fill((pim[0].momentum()+K0S[0].momentum()+K0L[0].momentum()).mass());
	    _h_3B_kspimkl[1]->fill((pim[0].momentum()+K0S[0].momentum()).mass());
	    _h_3B_kspimkl[2]->fill((K0S[0].momentum() +K0L[0].momentum()).mass());
	    _h_3B_kspimkl[3]->fill((pim[0].momentum()+K0L[0].momentum()).mass());
	  }
	}
	else if(nstable==5) {
	  if(pi0.size()==3&&pim.size()==1) {
	    _h_4B_pipi[0]     ->fill( (pi0[0].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipi[0]     ->fill( (pi0[1].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipi[0]     ->fill( (pi0[2].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipi[1]     ->fill( (pi0[0].momentum()+pi0[1].momentum()).mass());
	    _h_4B_pipi[1]     ->fill( (pi0[0].momentum()+pi0[2].momentum()).mass());
	    _h_4B_pipi[1]     ->fill( (pi0[1].momentum()+pi0[2].momentum()).mass());
	    _h_4B_pipipi[0]   ->fill( (pi0[0].momentum()+pi0[1].momentum()+pi0[2].momentum()).mass());
	    _h_4B_pipipi[1]   ->fill( (pi0[0].momentum()+pi0[1].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipipi[1]   ->fill( (pi0[0].momentum()+pi0[2].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipipi[1]   ->fill( (pi0[1].momentum()+pi0[2].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipipipi[0] ->fill( (pi0[0].momentum()+pi0[1].momentum()+pi0[2].momentum()+pim[0].momentum()).mass());
	  }
	  else if(pi0.size()==1&&pip.size()==1&&pim.size()==2) {
	    _h_4B_pipi[2] ->fill((pi0[0].momentum()+pip[0].momentum()).mass());
	    _h_4B_pipi[3] ->fill((pi0[0].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipi[3] ->fill((pi0[0].momentum()+pim[1].momentum()).mass());
	    _h_4B_pipi[4] ->fill((pip[0].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipi[4] ->fill((pip[0].momentum()+pim[1].momentum()).mass());
	    _h_4B_pipi[5] ->fill((pim[0].momentum()+pim[1].momentum()).mass());
	    _h_4B_pipipi[2]   ->fill( (pi0[0].momentum()+pip[0].momentum()+pim[0].momentum()).mass());
	    _h_4B_pipipi[2]   ->fill( (pi0[0].momentum()+pip[0].momentum()+pim[1].momentum()).mass());
	    _h_4B_pipipi[3]   ->fill( (pip[0].momentum()+pim[0].momentum()+pim[1].momentum()).mass());
	    _h_4B_pipipi[4]   ->fill( (pi0[0].momentum()+pim[0].momentum()+pim[1].momentum()).mass());
	    _h_4B_pipipipi[1] ->fill( (pi0[0].momentum()+pip[0].momentum()+pim[0].momentum()+pim[1].momentum()).mass());
	  }
	}
	else if(nstable==6) {
	  // 2 pi0 2pi- pi+
	  if(pi0.size()==2&&pim.size()==2&&pip.size()==1) {
	    FourMomentum ptotal = pim[0].momentum()+pim[1].momentum()+
	      pip[0].momentum()+pi0[0].momentum()+pi0[1].momentum();
	    _h_5B_pipi1[0]->fill((pim[0].momentum()+pim[1].momentum()).mass());
	    _h_5B_pipi1[1]->fill((pim[0].momentum()+pip[0].momentum()).mass());
	    _h_5B_pipi1[1]->fill((pim[1].momentum()+pip[0].momentum()).mass());
	    _h_5B_pipi1[2]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	    _h_5B_pipi1[2]->fill((pim[0].momentum()+pi0[1].momentum()).mass());
	    _h_5B_pipi1[2]->fill((pim[1].momentum()+pi0[0].momentum()).mass());
	    _h_5B_pipi1[2]->fill((pim[1].momentum()+pi0[1].momentum()).mass());
	    _h_5B_pipi1[3]->fill((pip[0].momentum()+pi0[0].momentum()).mass());
	    _h_5B_pipi1[3]->fill((pip[0].momentum()+pi0[1].momentum()).mass());
	    _h_5B_pipi1[4]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	    _h_5B_pipipi1[0]->fill((pim[0].momentum()+pim[1].momentum()-ptotal).mass());
	    _h_5B_pipipi1[1]->fill((pim[0].momentum()+pip[0].momentum()-ptotal).mass());
	    _h_5B_pipipi1[1]->fill((pim[1].momentum()+pip[0].momentum()-ptotal).mass());
	    _h_5B_pipipi1[2]->fill((pim[0].momentum()+pi0[0].momentum()-ptotal).mass());
	    _h_5B_pipipi1[2]->fill((pim[0].momentum()+pi0[1].momentum()-ptotal).mass());
	    _h_5B_pipipi1[2]->fill((pim[1].momentum()+pi0[0].momentum()-ptotal).mass());
	    _h_5B_pipipi1[2]->fill((pim[1].momentum()+pi0[1].momentum()-ptotal).mass());
	    _h_5B_pipipi1[3]->fill((pip[0].momentum()+pi0[0].momentum()-ptotal).mass());
	    _h_5B_pipipi1[3]->fill((pip[0].momentum()+pi0[1].momentum()-ptotal).mass());
	    _h_5B_pipipi1[4]->fill((pi0[0].momentum()+pi0[1].momentum()-ptotal).mass());	    
	    _h_5B_pipipipi1[0]->fill((ptotal-pim[0].momentum()).mass());
	    _h_5B_pipipipi1[0]->fill((ptotal-pim[1].momentum()).mass());
	    _h_5B_pipipipi1[1]->fill((ptotal-pip[0].momentum()).mass());
	    _h_5B_pipipipi1[2]->fill((ptotal-pi0[0].momentum()).mass());
	    _h_5B_pipipipi1[2]->fill((ptotal-pi0[1].momentum()).mass());
	    _h_5B_q1->fill(ptotal.mass());
	  }
	  // 4 pi0  pi-
	  else if(pi0.size()==4&&pim.size()==1) {
	    FourMomentum ptotal = pi0[0].momentum()+pi0[1].momentum()+pi0[2].momentum()+
	      pi0[3].momentum()+pim[0].momentum();
	    _h_5B_pipi2[0]->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	    _h_5B_pipi2[0]->fill((pim[0].momentum()+pi0[1].momentum()).mass());
	    _h_5B_pipi2[0]->fill((pim[0].momentum()+pi0[2].momentum()).mass());
	    _h_5B_pipi2[0]->fill((pim[0].momentum()+pi0[3].momentum()).mass());
	    _h_5B_pipi2[1]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	    _h_5B_pipi2[1]->fill((pi0[0].momentum()+pi0[2].momentum()).mass());
	    _h_5B_pipi2[1]->fill((pi0[0].momentum()+pi0[3].momentum()).mass());
	    _h_5B_pipi2[1]->fill((pi0[1].momentum()+pi0[2].momentum()).mass());
	    _h_5B_pipi2[1]->fill((pi0[1].momentum()+pi0[3].momentum()).mass());
	    _h_5B_pipi2[1]->fill((pi0[2].momentum()+pi0[3].momentum()).mass());
	    _h_5B_pipipi2[0]->fill((pim[0].momentum()+pi0[0].momentum()-ptotal).mass());
	    _h_5B_pipipi2[0]->fill((pim[0].momentum()+pi0[1].momentum()-ptotal).mass());
	    _h_5B_pipipi2[0]->fill((pim[0].momentum()+pi0[2].momentum()-ptotal).mass());
	    _h_5B_pipipi2[0]->fill((pim[0].momentum()+pi0[3].momentum()-ptotal).mass());
	    _h_5B_pipipi2[1]->fill((pi0[0].momentum()+pi0[1].momentum()-ptotal).mass());
	    _h_5B_pipipi2[1]->fill((pi0[0].momentum()+pi0[2].momentum()-ptotal).mass());
	    _h_5B_pipipi2[1]->fill((pi0[0].momentum()+pi0[3].momentum()-ptotal).mass());
	    _h_5B_pipipi2[1]->fill((pi0[1].momentum()+pi0[2].momentum()-ptotal).mass());
	    _h_5B_pipipi2[1]->fill((pi0[1].momentum()+pi0[3].momentum()-ptotal).mass());
	    _h_5B_pipipi2[1]->fill((pi0[2].momentum()+pi0[3].momentum()-ptotal).mass());
	    _h_5B_pipipipi2[0]->fill((ptotal-pim[0].momentum()).mass());
	    _h_5B_pipipipi2[1]->fill((ptotal-pi0[0].momentum()).mass());
	    _h_5B_pipipipi2[1]->fill((ptotal-pi0[1].momentum()).mass());
	    _h_5B_pipipipi2[1]->fill((ptotal-pi0[2].momentum()).mass());
	    _h_5B_pipipipi2[1]->fill((ptotal-pi0[3].momentum()).mass());
	    _h_5B_q2->fill(ptotal.mass());
	  }
	  // 3 pi- 2pi+
	  else if(pim.size()==3&&pip.size()==2) {
	    FourMomentum ptotal = pim[0].momentum()+pim[1].momentum()+
	      pim[2].momentum()+pip[0].momentum()+pip[1].momentum();
	    _h_5B_pipi3[0]->fill((pip[0].momentum()+pip[1].momentum()).mass());
	    _h_5B_pipi3[1]->fill((pim[0].momentum()+pip[0].momentum()).mass());
	    _h_5B_pipi3[1]->fill((pim[0].momentum()+pip[1].momentum()).mass());
	    _h_5B_pipi3[1]->fill((pim[1].momentum()+pip[0].momentum()).mass());
	    _h_5B_pipi3[1]->fill((pim[1].momentum()+pip[1].momentum()).mass());
	    _h_5B_pipi3[1]->fill((pim[2].momentum()+pip[0].momentum()).mass());
	    _h_5B_pipi3[1]->fill((pim[2].momentum()+pip[1].momentum()).mass());
	    _h_5B_pipi3[2]->fill((pim[0].momentum()+pim[1].momentum()).mass());
	    _h_5B_pipi3[2]->fill((pim[0].momentum()+pim[2].momentum()).mass());
	    _h_5B_pipi3[2]->fill((pim[1].momentum()+pim[2].momentum()).mass());
	    _h_5B_pipipi3[0]->fill((pip[0].momentum()+pip[1].momentum()-ptotal).mass());
	    _h_5B_pipipi3[1]->fill((pim[0].momentum()+pip[0].momentum()-ptotal).mass());
	    _h_5B_pipipi3[1]->fill((pim[0].momentum()+pip[1].momentum()-ptotal).mass());
	    _h_5B_pipipi3[1]->fill((pim[1].momentum()+pip[0].momentum()-ptotal).mass());
	    _h_5B_pipipi3[1]->fill((pim[1].momentum()+pip[1].momentum()-ptotal).mass());
	    _h_5B_pipipi3[1]->fill((pim[2].momentum()+pip[0].momentum()-ptotal).mass());
	    _h_5B_pipipi3[1]->fill((pim[2].momentum()+pip[1].momentum()-ptotal).mass());
	    _h_5B_pipipi3[2]->fill((pim[0].momentum()+pim[1].momentum()-ptotal).mass());
	    _h_5B_pipipi3[2]->fill((pim[0].momentum()+pim[2].momentum()-ptotal).mass());
	    _h_5B_pipipi3[2]->fill((pim[1].momentum()+pim[2].momentum()-ptotal).mass());
	    _h_5B_pipipipi3[0]->fill((ptotal-pim[0].momentum()).mass());
	    _h_5B_pipipipi3[0]->fill((ptotal-pim[1].momentum()).mass());
	    _h_5B_pipipipi3[0]->fill((ptotal-pim[2].momentum()).mass());
	    _h_5B_pipipipi3[1]->fill((ptotal-pip[0].momentum()).mass());
	    _h_5B_pipipipi3[1]->fill((ptotal-pip[1].momentum()).mass());
	    _h_5B_q3->fill(ptotal.mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // leptonic
      normalize(_h_2B_m2enu);
      normalize(_h_2B_menu );
      // 1 hadron
      normalize(_h_1B_xpi);
      // 2 hadrons
      normalize(_h_2B_m2pipi);
      normalize(_h_2B_mpipi );
      normalize(_h_2B_m2munu);
      normalize(_h_2B_mmunu );
      normalize(_h_2B_m2KpiA);
      normalize(_h_2B_mKpiA );
      normalize(_h_2B_m2KpiB);
      normalize(_h_2B_mKpiB );
      normalize(_h_2B_m2Keta);
      normalize(_h_2B_mKeta );
      normalize(_h_2B_m2KK  );
      normalize(_h_2B_mKK   );
      // 3 hadrons
      for(unsigned int ix=0;ix<4;++ix) {
	if(ix<3) {
	  normalize(_h_3B_pippimpim  [ix]);
	  normalize(_h_3B_pi0pi0pim  [ix]);
	  normalize(_h_3B_pi0pi0km   [ix]);
	  normalize(_h_3B_kspimks    [ix]);
	  normalize(_h_3B_klpimkl    [ix]);
	}
	normalize(_h_3B_kmpimkp    [ix]);
	normalize(_h_3B_kmpi0k0    [ix]);
	normalize(_h_3B_kmpimpip   [ix]);
	normalize(_h_3B_pimk0pi0   [ix]);
	normalize(_h_3B_pimpi0eta  [ix]);
	normalize(_h_3B_pimpi0gamma[ix]);
	normalize(_h_3B_kspimkl    [ix]);
      }
      // 4 pion decays
      for(unsigned int ix=0;ix<5;++ix) {
	normalize(_h_4B_pipi  [ix]);
	normalize(_h_4B_pipipi[ix]);
      }
      normalize(_h_4B_pipi[5]);
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_4B_pipipipi[ix]);
      }
      // 5 pions
      normalize(_h_5B_q1);
      for(unsigned int ix=0;ix<5;++ix) {
	normalize(_h_5B_pipi1);
	normalize(_h_5B_pipipi1);
      }
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h_5B_pipipipi1);
      }
      // 4 pi0  pi-
      normalize(_h_5B_q2);
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_5B_pipi2);
	normalize(_h_5B_pipipi2);
	normalize(_h_5B_pipipipi2);
      }
      // 3 pi- 2 pi+
      normalize(_h_5B_q3);
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h_5B_pipi3);
	normalize(_h_5B_pipipi3);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	normalize(_h_5B_pipipipi3);
      }
    }

    //@}


    /// @name Histograms
    //@{
    // histograms for leptonic decay
    Histo1DPtr _h_2B_m2enu,_h_2B_menu;
    Histo1DPtr _h_2B_m2munu,_h_2B_mmunu;
    // histograms for 1 hadron decay
    Histo1DPtr _h_1B_xpi;
    // histograms for 2 hadron decay
    Histo1DPtr _h_2B_m2pipi,_h_2B_mpipi;
    Histo1DPtr _h_2B_m2KpiA,_h_2B_m2KpiB,_h_2B_mKpiA,_h_2B_mKpiB;
    Histo1DPtr _h_2B_m2Keta,_h_2B_mKeta;
    Histo1DPtr _h_2B_m2KK,_h_2B_mKK;
    // histograms for 3 hadronc decay
    //  Histograms for tau^- -> nu_tau pi^+pi^-pi^-  
    vector<Histo1DPtr> _h_3B_pippimpim;
    //Histograms for tau^- -> nu_tau pi^0pi^0pi^-  
    vector<Histo1DPtr> _h_3B_pi0pi0pim;
    //  Histograms for tau^- -> nu_tau K^-K^+pi^-      
    vector<Histo1DPtr> _h_3B_kmpimkp;
    //  Histograms for tau^- -> nu_tau K^-K^0pi^0      
    vector<Histo1DPtr> _h_3B_kmpi0k0;
    //  Histograms for tau^- -> nu_tau pi^0pi^0K^-   
    vector<Histo1DPtr> _h_3B_pi0pi0km;
    //  Histograms for tau^- -> nu_tau K^-pi^-pi^+    
    vector<Histo1DPtr> _h_3B_kmpimpip; 
    //  Histograms for tau^- -> nu_tau pi^-K^0pi^0    
    vector<Histo1DPtr> _h_3B_pimk0pi0;
    //  Histograms for tau^- -> nu_tau pi^-pi^0eta   
    vector<Histo1DPtr> _h_3B_pimpi0eta;
    //  Histograms for tau^- -> nu_tau pi^-pi^0gamma 
    vector<Histo1DPtr> _h_3B_pimpi0gamma;
    //  Histograms for tau^- -> nu_tau K^0_SK^0_Spi^-
    vector<Histo1DPtr> _h_3B_kspimks;
    //  Histograms for tau^- -> nu_tau K^0_LK^0_Lpi^-
    vector<Histo1DPtr> _h_3B_klpimkl;
    //  Histograms for tau^- -> nu_tau K^0_SK^0_Lpi^-
    vector<Histo1DPtr> _h_3B_kspimkl;
    // histograms for 4 pion decay
    //  Histograms for the pipi mass distributions
    vector<Histo1DPtr>  _h_4B_pipi;
    //  Histograms for the pipipi mass distributions
    vector<Histo1DPtr>  _h_4B_pipipi;
    //  Histograms for the pipipipi mass distributions
    vector<Histo1DPtr>  _h_4B_pipipipi;
    // histograms for 5 pion decay
    // 2 pi0 2 pi- pi+
    Histo1DPtr _h_5B_q1;
    vector<Histo1DPtr> _h_5B_pipi1;
    vector<Histo1DPtr> _h_5B_pipipi1;
    vector<Histo1DPtr> _h_5B_pipipipi1;
    // 4 pi0 pi-
    Histo1DPtr _h_5B_q2;
    vector<Histo1DPtr> _h_5B_pipi2;
    vector<Histo1DPtr> _h_5B_pipipi2;
    vector<Histo1DPtr> _h_5B_pipipipi2;
    // 3 pi- 2 pi+
    Histo1DPtr _h_5B_q3;
    vector<Histo1DPtr> _h_5B_pipi3;
    vector<Histo1DPtr> _h_5B_pipipi3;
    vector<Histo1DPtr> _h_5B_pipipipi3;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TAU_Decay);


}
