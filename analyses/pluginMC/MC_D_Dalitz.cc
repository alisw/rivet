// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Monte Carlo analysis of D meson dalitz decays
  class MC_D_Dalitz : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_D_Dalitz);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411 or
						Cuts::abspid==421 or
						Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");

      // Book histograms
      book(_h_plus1, "h_plus1"   ,200,0.,3. );
      book(_h_minus1, "h_minus1"  ,200,0.,3.2 );
      book(_h_pipi1, "h_pipi1"   ,200,0.,2. );
      book(_h_minus2, "h_minus2"  ,200,0.,3.2 );
      book(_h_neutral2, "h_neutral2",200,0.,3.2 );
      book(_h_pipi2, "h_pipi2"   ,200,0.,2. );
      book(_h_Kpilow3, "h_Kpilow3"  ,200,0.,2. );
      book(_h_Kpihigh3, "h_Kpihigh3" ,200,0.,3.2 );
      book(_h_Kpiall3, "h_Kpiall3"  ,200,0.,3. );
      book(_h_pipi3, "h_pipi3"    ,200,0.,2. );
      book(_h_Kpip4, "h_Kpip4"   ,200,0.,3.2 );
      book(_h_pipi4, "h_pipi4"   ,200,0.,2. );
      book(_h_Kpi04, "h_Kpi04"   ,200,0.,3.2);
      book(_h_kppim5, "h_kppim5"   ,200,0.,3. );
      book(_h_kppip5, "h_kppip5"   ,200,0.,3.1 );
      book(_h_pippim5, "h_pippim5"  ,200,0.,2. );
      book(_h_kppim6, "h_kppim6"   ,200,0.,3.5);
      book(_h_kppip6, "h_kppip6"   ,200,0.,3.5);
      book(_h_pippim6, "h_pippim6"  ,200,0.,2.5);
      book(_h_kpkm1 ,"h_kpkm1",200,0.9,3.5);
      book(_h_kppip7,"h_kppip7",200,0.3,3.5);
      book(_h_kmpip1,"h_kmpip1",200,0.3,3.5);
      book(_h_pipi5, "h_pipi5"   ,200,0.,3.2 );
      book(_h_pipi6, "h_pipi6"   ,200,0.,3.2 );
      book(_h_pipi7, "h_pipi7"   ,200,0.,3.2 );
      book(_dalitz1, "dalitz1"    ,50,0.3,3.2,50,0.3,3.2);
      book(_dalitz2, "dalitz2"    ,50,0.3,3. ,50,0.3,3. );
      book(_dalitz3, "dalitz3"    ,50,0.3,2. ,50,0.07,2. );
      book(_dalitz4, "dalitz4"    ,50,0.3,3.1 ,50,0.07,2. );
      book(_dalitz5, "dalitz5"    ,50,0.,3. ,50,0.,2. );
      book(_dalitz6, "dalitz6"    ,50,0.3,3.5,50,0.07,2.5);
      book(_dalitz7, "dalitz7"    ,50,0.3,3.5,50,0.07,2.5);
      book(_dalitz8, "dalitz8"    ,50,0.,3.2,50,0.,3.2);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 211,1}, {-211,1}, { 310,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 211,1}, {-321,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode2CC = { {-211,1}, { 321,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode3   = { { 211,1}, {-211,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode4   = { { 211,2}, {-321,1} };
      static const map<PdgId,unsigned int> & mode4CC = { {-211,2}, { 321,1} };
      static const map<PdgId,unsigned int> & mode5   = { { 211,1}, { 111,1}, { 310,1} };
      static const map<PdgId,unsigned int> & mode5CC = { {-211,1}, { 111,1}, { 310,1} };
      static const map<PdgId,unsigned int> & mode6   = { { 211,1}, {-211,1}, { 321,1} };
      static const map<PdgId,unsigned int> & mode6CC = { { 211,1}, {-211,1}, {-321,1} };
      static const map<PdgId,unsigned int> & mode7   = { { 321,1}, {-321,1}, { 211,1} };
      static const map<PdgId,unsigned int> & mode7CC = { { 321,1}, {-321,1}, {-211,1} };
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      for(unsigned int ix=0;ix<DD.decaying().size();++ix) {
	int sign = DD.decaying()[ix].pid()/DD.decaying()[ix].abspid();
	if(DD.decaying()[ix].abspid()==421) {
	  if ( DD.modeMatches(ix,3,mode1  )) {
	    const Particle & K0  = DD.decayProducts()[ix].at(      310)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	    double mminus = (pim.momentum()+K0.momentum() ).mass2();
	    double mplus  = (pip.momentum()+K0.momentum() ).mass2();
	    double mpipi  = (pip.momentum()+pim.momentum()).mass2();
	    _h_plus1  ->fill(mplus);
	    _h_minus1 ->fill(mminus);
	    _h_pipi1  ->fill(mpipi);
	    _dalitz1   ->fill(mplus,mminus);
	  }
	  else if( ( DD.decaying()[ix].pid()>0 && DD.modeMatches(ix,3,mode2  )) ||
		   ( DD.decaying()[ix].pid()<0 && DD.modeMatches(ix,3,mode2CC))) {
	    const Particle & pi0 = DD.decayProducts()[ix].at(      111)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & Km  = DD.decayProducts()[ix].at(-sign*321)[0];
	    double mneut  = (Km.momentum()+pip.momentum()).mass2();
	    double mminus = (Km.momentum()+pi0.momentum()).mass2();
	    double mpipi  = (pip.momentum()+pi0.momentum()).mass2();
	    _h_neutral2  ->fill(mneut);
	    _h_minus2 ->fill(mminus);
	    _h_pipi2  ->fill(mpipi);
	    _dalitz2->fill(mminus,mneut);
	  }
	  else if(DD.modeMatches(ix,3,mode3)) {
	    const Particle & pi0 = DD.decayProducts()[ix].at(      111)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	    double mneut  = (pim.momentum()+pip.momentum()).mass2();
	    double mminus = (pim.momentum()+pi0.momentum()).mass2();
	    double mplus  = (pip.momentum()+pi0.momentum()).mass2();
	    _h_pipi5 ->fill(mplus);
	    _h_pipi6 ->fill(mminus);
	    _h_pipi7 ->fill(mneut);
	    _dalitz8 ->fill(mplus,mminus);
	  }
	}
	else if(DD.decaying()[ix].abspid()==411) {
	  if(DD.modeMatches(ix,3,mode4  ) || DD.modeMatches(ix,3,mode4CC)) {
	    const Particles & pip = DD.decayProducts()[ix].at( sign*211);
	    const Particle  & Km  = DD.decayProducts()[ix].at(-sign*321)[0];
	    double mplus  = (Km.momentum() +pip[0].momentum()).mass2();
	    double mminus = (Km.momentum() +pip[1].momentum()).mass2();
	    double mpipi  = (pip[0].momentum()+pip[1].momentum()).mass2();
	    if(mplus<mminus) swap(mplus,mminus);
	    _h_Kpilow3 ->fill( mminus);
	    _h_Kpihigh3->fill( mplus);
	    _h_Kpiall3 ->fill( mminus);
	    _h_Kpiall3 ->fill( mplus);
	    _h_pipi3   ->fill( mpipi);
	    _dalitz3->fill(mminus,mpipi);
	  }
	  else if(DD.modeMatches(ix,3,mode5  ) || DD.modeMatches(ix,3,mode5CC)) {
	    const Particle & pi0 = DD.decayProducts()[ix].at(      111)[0];
	    const Particle & K0  = DD.decayProducts()[ix].at(      310)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    double mminus = (K0.momentum()+pip.momentum()).mass2();
	    double mplus  = (K0.momentum()+pi0.momentum()).mass2();
	    double mpipi  = (pip.momentum()+pi0.momentum()).mass2();
	    _h_Kpip4 ->fill( mminus);
	    _h_pipi4 ->fill( mpipi );
	    _h_Kpi04 ->fill( mplus );
	    _dalitz4->fill(mplus,mpipi);
	  }
	  else if(DD.modeMatches(ix,3,mode6  ) || DD.modeMatches(ix,3,mode6CC)) {
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	    const Particle & Kp  = DD.decayProducts()[ix].at( sign*321)[0];
	    double mplus  = (Kp .momentum()+pip.momentum()).mass2();
	    double mminus = (Kp .momentum()+pim.momentum()).mass2();
	    double mpipi  = (pip.momentum()+pim.momentum()).mass2();
	    _h_kppim5 ->fill(mminus);
	    _h_kppip5 ->fill(mplus );
	    _h_pippim5->fill(mpipi );
	    _dalitz5->fill(mminus,mpipi);
	  }
	}
	else if(DD.decaying()[ix].abspid()==431) {
	  if(DD.modeMatches(ix,3,mode6  ) || DD.modeMatches(ix,3,mode6CC)) {
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    const Particle & pim = DD.decayProducts()[ix].at(-sign*211)[0];
	    const Particle & Kp  = DD.decayProducts()[ix].at( sign*321)[0];
	    double mplus  = (Kp .momentum()+pip.momentum()).mass2();
	    double mminus = (Kp .momentum()+pim.momentum()).mass2();
	    double mpipi  = (pip.momentum()+pim.momentum()).mass2();
	    _h_kppim6 ->fill(mminus);
	    _h_kppip6 ->fill(mplus );
	    _h_pippim6->fill(mpipi );
	    _dalitz6->fill(mminus,mpipi);
	  }
	  else if(DD.modeMatches(ix,3,mode7  ) || DD.modeMatches(ix,3,mode7CC)) {
	    const Particle & Kp  = DD.decayProducts()[ix].at( sign*321)[0];
	    const Particle & Km  = DD.decayProducts()[ix].at(-sign*321)[0];
	    const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
	    double mplus  = (Kp.momentum()+pip.momentum()).mass2();
	    double mminus = (Km.momentum()+pip.momentum()).mass2();
	    double mKK    = (Kp.momentum()+Km .momentum()).mass2();
	    _h_kpkm1->fill(mKK);
	    _h_kppip7->fill(mplus);
	    _h_kmpip1->fill(mminus);
	    _dalitz7->fill(mKK,mminus);
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_plus1);
      normalize(_h_minus1);
      normalize(_h_pipi1);
      normalize(_dalitz1);
      normalize(_h_minus2);
      normalize(_h_pipi2);
      normalize(_h_neutral2);
      normalize(_dalitz2);
      normalize(_h_Kpilow3);
      normalize(_h_Kpihigh3);
      normalize(_h_Kpiall3);
      normalize(_h_pipi3);
      normalize(_dalitz3);
      normalize(_h_Kpip4);
      normalize(_h_pipi4);
      normalize(_h_Kpi04);
      normalize(_dalitz4);
      normalize(_h_kppim5);
      normalize(_h_kppip5);
      normalize(_h_pippim5);
      normalize(_dalitz5);
      normalize(_h_kppim6);
      normalize(_h_kppip6);
      normalize(_h_pippim6);
      normalize(_dalitz6);
      normalize(_h_kpkm1);
      normalize(_h_kppip7);
      normalize(_h_kmpip1);
      normalize(_dalitz7);
      normalize(_h_pipi5);
      normalize(_h_pipi6);
      normalize(_h_pipi7);
      normalize(_dalitz8);
    }
    //@}

    /// @name Histograms
    //@{
    // Histograms for D^0\to \bar{K}^0\pi^+\pi^-
    //m^2_+
    Histo1DPtr _h_plus1;
    //m^2_+
    Histo1DPtr _h_minus1;
    //m^2_{\pi\pi}
    Histo1DPtr _h_pipi1;
    // Dalitz plot
    Histo2DPtr _dalitz1;
    
    // Histograms for D^0\to K^-\pi^+\pi^0
    // Histogram for the K^-\pi^+ mass
    Histo1DPtr _h_minus2;
    // Histogram for the \pi^+\pi^0 mass
    Histo1DPtr _h_pipi2;
    // Histogram for the K^-\pi^0 mass
    Histo1DPtr _h_neutral2;
    // Dalitz plot
    Histo2DPtr _dalitz2;

    // Histograms for D^+\to K^-\pi^+\pi^+
    // Histogram for K^-\pi^+ low
    Histo1DPtr _h_Kpilow3;
    // Histogram for K^-\pi^+ high
    Histo1DPtr _h_Kpihigh3;
    // Histogram for K^-\pi^+ all
    Histo1DPtr _h_Kpiall3;
    // Histogram for \pi^+\pi^-
    Histo1DPtr _h_pipi3;
    // Dalitz plot
    Histo2DPtr _dalitz3;

    // Histograms for D^+\to\bar{K}^0\pi^+\pi^0
    // Histogram for the \bar{K}^0\pi^+ mass
    Histo1DPtr _h_Kpip4;
    // Histogram for the \pi^+\pi^0 mass
    Histo1DPtr _h_pipi4;
    // Histogram for the \bar{K}^0\pi^0 mass
    Histo1DPtr _h_Kpi04;
    // Dalitz plot
    Histo2DPtr _dalitz4;

    // Histograms for D^+\to K^+\pi^-\pi^+
    // Histogram for K^+\pi^-
    Histo1DPtr _h_kppim5;
    // Histogram for K^+\pi^+
    Histo1DPtr _h_kppip5;
    // Histogram for \pi^+\pi^-
    Histo1DPtr _h_pippim5;
    // Dalitz plot
    Histo2DPtr _dalitz5;

    // Histograms for D_s^+\to K^+\pi^-\pi^+
    // Histogram for K^+\pi^-
    Histo1DPtr _h_kppim6;
    // Histogram for K^+\pi^+
    Histo1DPtr _h_kppip6;
    // Histogram for \pi^+\pi^-
    Histo1DPtr _h_pippim6;
    // Dalitz plot
    Histo2DPtr _dalitz6;

    // Histograms for D_s^+\to K^+K^-\pi^+
    // Histogram for K^+K^-
    Histo1DPtr _h_kpkm1;
    // Histogram for K^+\pi^+
    Histo1DPtr _h_kppip7;
    // Histogram for K^-\pi^+
    Histo1DPtr _h_kmpip1;
    // Dalitz plot
    Histo2DPtr _dalitz7;

    // Histograms for D0 -> pi+pi-pi0
    // Histogram for pi+pi0
    Histo1DPtr _h_pipi5;
    // Histogram for pi-pi0
    Histo1DPtr _h_pipi6;
    // Histogram for pi+pi-
    Histo1DPtr _h_pipi7;
    // Dalitz plot
    Histo2DPtr _dalitz8;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_D_Dalitz);


}
