// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta_c Dalitz decays
  class BABAR_2021_I1867843 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2021_I1867843);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid== 441);
      declare(ufs, "UFS");
      DecayedParticles ETAC(ufs);
      ETAC.addStable(PID::PI0);
      ETAC.addStable(PID::K0S);
      ETAC.addStable(PID::ETA);
      ETAC.addStable(PID::ETAPRIME);
      declare(ETAC,"ETAC");
      // histograms
      book(_h_KK    ,1,1,1);
      book(_h_etaPK ,1,1,2);
      book(_h_pipi1 ,2,1,1);
      book(_h_etaPpi,2,1,2);
      book(_h_pipi2 ,3,1,1);
      book(_h_etapi ,3,1,2);
      book(_dalitz1, "dalitz1",50,2.,6.5,50,2.,6.5);
      book(_dalitz2, "dalitz2",50,0.,9. ,50,0.,9. );
      book(_dalitz3, "dalitz3",50,0.,8. ,50,0.,8. );
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1  = { { 321,1}, {-321,1}, { 331,1}};
      static const map<PdgId,unsigned int> & mode2  = { { 211,1}, {-211,1}, { 331,1}};
      static const map<PdgId,unsigned int> & mode3  = { { 211,1}, {-211,1}, { 221,1}};
      DecayedParticles ETAC = apply<DecayedParticles>(event, "ETAC");
      // loop over particles
      for(unsigned int ix=0;ix<ETAC.decaying().size();++ix) {
	//K+ K- eta'
	if (ETAC.modeMatches(ix,3,mode1)&&
	    ETAC.decaying()[ix].mass()>2.93 && ETAC.decaying()[ix].mass()<3.03) {
	  const Particle & Kp   = ETAC.decayProducts()[ix].at( 321)[0];
	  const Particle & Km   = ETAC.decayProducts()[ix].at(-321)[0];
	  const Particle & etaP = ETAC.decayProducts()[ix].at( 331)[0];
	  double mplus  = (Kp.momentum()+etaP.momentum()).mass2();
	  double mminus = (Km.momentum()+etaP.momentum()).mass2();
	  double mKK    = (Kp.momentum()+Km .momentum()).mass2();
	  _h_KK   ->fill(sqrt(mKK));
	  _h_etaPK->fill(sqrt(mplus));
	  _h_etaPK->fill(sqrt(mminus));
	  _dalitz1->fill(mplus,mminus);
	}
	// pi+ pi- eta'
	else if (ETAC.modeMatches(ix,3,mode2)&&
		 ETAC.decaying()[ix].mass()>2.93 && ETAC.decaying()[ix].mass()<3.03) {
	  const Particle & pip  = ETAC.decayProducts()[ix].at( 211)[0];
	  const Particle & pim  = ETAC.decayProducts()[ix].at(-211)[0];
	  const Particle & etaP = ETAC.decayProducts()[ix].at( 331)[0];
	  double mplus  = (pip.momentum()+etaP.momentum()).mass2();
	  double mminus = (pim.momentum()+etaP.momentum()).mass2();
	  double mpipi    = (pip.momentum()+pim .momentum()).mass2();
	  _h_pipi1 ->fill(sqrt(mpipi));
	  _h_etaPpi->fill(sqrt(mplus));
	  _h_etaPpi->fill(sqrt(mminus));
	  _dalitz2->fill(mplus,mminus);
	}
	// pi+ pi- eta
	else if (ETAC.modeMatches(ix,3,mode3)&&
		 ETAC.decaying()[ix].mass()>2.92 && ETAC.decaying()[ix].mass()<3.02) {
	  const Particle & pip = ETAC.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = ETAC.decayProducts()[ix].at(-211)[0];
	  const Particle & eta = ETAC.decayProducts()[ix].at( 221)[0];
	  double mplus  = (pip.momentum()+eta.momentum()).mass2();
	  double mminus = (pim.momentum()+eta.momentum()).mass2();
	  double mpipi    = (pip.momentum()+pim .momentum()).mass2();
	  _h_pipi2->fill(sqrt(mpipi));
	  _h_etapi->fill(sqrt(mplus));
	  _h_etapi->fill(sqrt(mminus));
	  _dalitz3->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_KK    );
      normalize(_h_etaPK );
      normalize(_h_pipi1 );
      normalize(_h_etaPpi);
      normalize(_h_pipi2 );
      normalize(_h_etapi );
      normalize(_dalitz1 );
      normalize(_dalitz2 );
      normalize(_dalitz3 );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_KK,_h_etaPK;
    Histo1DPtr _h_pipi1,_h_etaPpi;
    Histo1DPtr _h_pipi2,_h_etapi;
    Histo2DPtr _dalitz1,_dalitz2,_dalitz3;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2021_I1867843);

}
