// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MC_OmegaPhia1_3Pion_Decay : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_OmegaPhia1_3Pion_Decay);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // book histograms a_1
      // Histograms for a_10 -> pi0pi0pi0
      book(_hist0, "hist0",200,0.2,1.5);
      // // dalitz plot
      book(_dalitz0, "dalitz0",50,0.2,1.5,50,0.2,1.5);
      // Histograms for a_1+ -> pi0pi0pi+
      // Mass of the pi0pi0 pair 
      book(_hist1A, "hist1A",200,0.2,1.5);
      // Mass of the pi0pi+ pair
      book(_hist1B, "hist1B",200,0.2,1.5);
      // dalitz plot
      book(_dalitz1, "dalitz1",50,0.2,1.5,50,0.2,1.5);
      // Histograms for a_10 -> pi+pi-pi0
      // Mass of the pi+pi- pair 
      book(_hist2A, "hist2A",200,0.2,1.5);
      // Mass of the pi+pi0 pair
      book(_hist2B, "hist2B",200,0.2,1.5);
      // Mass of the pi-pi0 pair
      book(_hist2C, "hist2C",200,0.2,1.5);
      // dalitz plot
      book(_dalitz2, "dalitz2",50,0.2,1.5,50,0.2,1.5);
      //  Histograms for a_1+ -> pi+pi+pi-
      // Mass of the pi+pi+ pair
      book(_hist3A, "hist3A",200,0.2,1.5);
      // Mass of the pi+pi- pair
      book(_hist3B, "hist3B",200,0.2,1.5);
      // dalitz plot
      book(_dalitz3, "dalitz3",50,0.2,1.5,50,0.2,1.5);
      
      // Book histograms omega/phi
      for(unsigned int ix=0;ix<2;++ix) {
	double mmax = ix==0 ? 0.8 : 1.0; 
	std::ostringstream title1; title1 << "xhist_" << ix+1;
	_h_xhist  .push_back(Histo1DPtr());
        book(_h_xhist.back(), title1.str(),200,-300.,300. );
	std::ostringstream title2; title2 << "yhist_" << ix+1;
	_h_yhist  .push_back(Histo1DPtr());
        book(_h_yhist.back(), title2.str(),200,0.   ,400. );
	std::ostringstream title3; title3 << "mplus_" << ix+1;
	_h_mplus  .push_back(Histo1DPtr());
        book(_h_mplus.back(), title3.str(),200,200.,mmax*1000.);
	std::ostringstream title4; title4 << "mminus_" << ix+1;
	_h_mminus .push_back(Histo1DPtr());
        book(_h_mminus.back(), title4.str(),200,200.,mmax*1000.);
	std::ostringstream title5; title5 << "m0_" << ix+1;
	_h_m0     .push_back(Histo1DPtr());
        book(_h_m0.back(), title5.str(),200,200.,mmax*1000.);
	std::ostringstream title6; title6 << "dalitz_" << ix+1;
	_h_dalitz.push_back(Histo2DPtr());
        book(_h_dalitz.back(), title6.str(),50,0.2,mmax,50,0.2,mmax);
      }
    }


    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles & pip , Particles & pim , Particles & pi0) {
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
	else if (id == PID::PI0) {
	  pi0.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p, nstable, pip, pim, pi0);
	}
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& meson :
	    apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::PHI || Cuts::pid==PID::OMEGA ||
							     Cuts::abspid==20213 || Cuts::pid==20113 )) {
      	unsigned int nstable(0);
       	Particles pip, pim, pi0;
      	findDecayProducts(meson, nstable, pip, pim, pi0);
	if(nstable !=3) continue;
	if(meson.pid()<0) {
	  swap(pim,pip);
	}
	if(meson.pid()== PID::PHI || meson.pid()==PID::OMEGA) {
	  if(pip.size()!=1 || pim.size()!=1 || pi0.size()!=1) continue;
	  unsigned int iloc = meson.pid() == PID::OMEGA ? 0 : 1;
	  LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(meson.momentum().betaVec());
	  FourMomentum pp = boost.transform(pip[0].momentum());
	  FourMomentum pm = boost.transform(pim[0].momentum());
	  FourMomentum p0 = boost.transform(pi0[0].momentum());
	  double mp = (pp+p0).mass(), mm  = (pm+pp).mass();
	  _h_mplus [iloc]->fill(mp/MeV);
	  _h_mminus[iloc]->fill((pm+p0).mass()/MeV);
	  _h_m0    [iloc]->fill(mm/MeV);
	  double x = pp.t()-pm.t();
	  double y = p0.t()-p0.mass();
	  _h_xhist[iloc]->fill(x/MeV);
	  _h_yhist[iloc]->fill(y/MeV);
	  _h_dalitz[iloc]->fill(mp,mm);
	}
	else {
	  // a_1+ -> pi+pi+pi-
	  if(pip.size()==2&&pim.size()==1) {
	    _hist3A->fill((pip[0].momentum()+pip[1].momentum()).mass());
	    _hist3B->fill((pip[0].momentum()+pim[0].momentum()).mass());
	    _hist3B->fill((pip[1].momentum()+pim[0].momentum()).mass());
	    _dalitz3->fill((pip[0].momentum()+pim[0].momentum()).mass(),(pip[1].momentum()+pim[0].momentum()).mass());
	    _dalitz3->fill((pip[1].momentum()+pim[0].momentum()).mass(),(pip[0].momentum()+pim[0].momentum()).mass());
	  }
	  // a_1+ -> pi0pi0pi+
	  else if(pip.size()==1&&pi0.size()==2) {
	    _hist1A->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	    _hist1B->fill((pip[0].momentum()+pi0[0].momentum()).mass());
	    _hist1B->fill((pip[0].momentum()+pi0[1].momentum()).mass());
	    _dalitz1->fill((pip[0].momentum()+pi0[0].momentum()).mass(),(pip[0].momentum()+pi0[1].momentum()).mass());
	    _dalitz1->fill((pip[0].momentum()+pi0[1].momentum()).mass(),(pip[0].momentum()+pi0[0].momentum()).mass());
	  }
	  // a_10 -> pi0pi0pi0
	  else if(pi0.size()==3) { 
	    _hist0->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	    _hist0->fill((pi0[0].momentum()+pi0[2].momentum()).mass());
	    _hist0->fill((pi0[1].momentum()+pi0[2].momentum()).mass());
	    _dalitz0->fill((pi0[0].momentum()+pi0[1].momentum()).mass(),(pi0[0].momentum()+pi0[2].momentum()).mass());
	    _dalitz0->fill((pi0[0].momentum()+pi0[1].momentum()).mass(),(pi0[1].momentum()+pi0[2].momentum()).mass());
	    _dalitz0->fill((pi0[0].momentum()+pi0[2].momentum()).mass(),(pi0[1].momentum()+pi0[2].momentum()).mass());
	    _dalitz0->fill((pi0[0].momentum()+pi0[2].momentum()).mass(),(pi0[0].momentum()+pi0[1].momentum()).mass());
	    _dalitz0->fill((pi0[1].momentum()+pi0[2].momentum()).mass(),(pi0[0].momentum()+pi0[1].momentum()).mass());
	    _dalitz0->fill((pi0[1].momentum()+pi0[2].momentum()).mass(),(pi0[0].momentum()+pi0[2].momentum()).mass());
	  }
	  // a_10 -> pi+pi-pi0
	  else if(pi0.size()==1&&pip.size()==1&&pim.size()==1) {
	    _hist2A->fill((pim[0].momentum()+pip[0].momentum()).mass());
	    _hist2B->fill((pip[0].momentum()+pi0[0].momentum()).mass());
	    _hist2C->fill((pim[0].momentum()+pi0[0].momentum()).mass());
	    _dalitz2->fill((pim[0].momentum()+pi0[0].momentum()).mass(),(pip[0].momentum()+pi0[0].momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // a_1
      normalize(_hist0);
      normalize(_dalitz0);
      normalize(_hist1A);
      normalize(_hist1B);
      normalize(_dalitz1);
      normalize(_hist2A);
      normalize(_hist2B);
      normalize(_hist2C);
      normalize(_dalitz2);
      normalize(_hist3A);
      normalize(_hist3B);
      normalize(_dalitz3);
      // omega/phi
      for(unsigned int ix=0;ix<2;++ix) {
        normalize(_h_xhist [ix]);
        normalize(_h_yhist [ix]);
        normalize(_h_mplus [ix]);
        normalize(_h_mminus[ix]);
        normalize(_h_m0    [ix]);
        normalize(_h_dalitz[ix]);
      }
    }

    //@}

    /// @name Histograms a_1
    //@{
    // Histograms for a_10 -> pi0pi0pi0
    Histo1DPtr _hist0;
    // dalitz plot
    Histo2DPtr _dalitz0;
    // Histograms for a_1+ -> pi0pi0pi+
    // Mass of the pi0pi0 pair 
    Histo1DPtr _hist1A;
    // Mass of the pi0pi+ pair
    Histo1DPtr _hist1B;
    // dalitz plot
    Histo2DPtr _dalitz1;
    // Histograms for a_10 -> pi+pi-pi0
    // Mass of the pi+pi- pair 
    Histo1DPtr _hist2A;
    // Mass of the pi+pi0 pair
    Histo1DPtr _hist2B;
    // Mass of the pi-pi0 pair
    Histo1DPtr _hist2C;
    // dalitz plot
    Histo2DPtr _dalitz2;
    //  Histograms for a_1+ -> pi+pi+pi-
    // Mass of the pi+pi+ pair
    Histo1DPtr _hist3A;
    // Mass of the pi+pi- pair
    Histo1DPtr _hist3B;
    // dalitz plot
    Histo2DPtr _dalitz3;
    //@}

    /// @name Histograms omega/phi
    //@{
    // Histogram for the x-values
    vector<Histo1DPtr> _h_xhist;
    // Histogram for the y-values
    vector<Histo1DPtr> _h_yhist;
    //  The mass of the \rho^+
    vector<Histo1DPtr> _h_mplus;
    //  The mass of the \rho^-
    vector<Histo1DPtr> _h_mminus;
    // The mass of the \rho^0
    vector<Histo1DPtr> _h_m0;
    // Dalitz plot
    vector<Histo2DPtr> _h_dalitz;
    //@}
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_OmegaPhia1_3Pion_Decay);


}
