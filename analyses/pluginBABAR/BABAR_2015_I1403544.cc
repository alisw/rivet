// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta_c -> K+K-pi0 KS0 K+-pi-+
  class BABAR_2015_I1403544 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2015_I1403544);


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
      book(_h_Kppip,3,1,1);
      book(_h_K0pip,3,1,2);
      book(_h_K0Kp ,3,1,3);
      book(_dalitz[0], "dalitz_1",50,0.,8.,50,0.,8.);
      book(_h_Kppi0,4,1,1);
      book(_h_Kmpi0,4,1,2);
      book(_h_KpKm ,4,1,3);
      book(_dalitz[1], "dalitz_2",50,0.,8.,50,0.,8.);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1}, {-321,1}, { 111,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 321,1}, {-211,1}, { 310,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-321,1}, { 211,1}, { 310,1}};
      DecayedParticles ETAC = apply<DecayedParticles>(event, "ETAC");
      // loop over particles
      for(unsigned int ix=0;ix<ETAC.decaying().size();++ix) {
	// K+ K- pi0
	if (ETAC.modeMatches(ix,3,mode1)&&
	    ETAC.decaying()[ix].mass()>2.922 && ETAC.decaying()[ix].mass()<3.036) {
	  const Particle & Kp  = ETAC.decayProducts()[ix].at( 321)[0];
	  const Particle & Km  = ETAC.decayProducts()[ix].at(-321)[0];
	  const Particle & pi0 = ETAC.decayProducts()[ix].at( 111)[0];
	  double mplus  = (Kp.momentum()+pi0.momentum()).mass2();
	  double mminus = (Km.momentum()+pi0.momentum()).mass2();
	  double mKK    = (Kp.momentum()+Km .momentum()).mass2();
	  _h_KpKm->fill(mKK);
	  _h_Kppi0->fill(mplus);
	  _h_Kmpi0->fill(mminus);
	  _dalitz[1]->fill(mplus,mminus);
	  continue;
	}
	else if(ETAC.decaying()[ix].mass()>2.922 && ETAC.decaying()[ix].mass()<3.039) {
	  int sign=1;
	  if     (ETAC.modeMatches(ix,3,mode2  )) sign= 1;
	  else if(ETAC.modeMatches(ix,3,mode2CC)) sign=-1;
	  else continue;
	  const Particle & KS0 = ETAC.decayProducts()[ix].at( 310)[0];
	  const Particle & Kp  = ETAC.decayProducts()[ix].at( sign*321)[0];
	  const Particle & pim = ETAC.decayProducts()[ix].at(-sign*211)[0];
	  double mplus  = (Kp.momentum()  + pim.momentum()).mass2();
	  double mminus = (KS0.momentum() + pim.momentum()).mass2();
	  double mKK    = (Kp.momentum()  + KS0.momentum()).mass2();
	  _h_K0Kp ->fill(mKK);
	  _h_Kppip->fill(mplus);
	  _h_K0pip->fill(mminus);
	  _dalitz[0]->fill(mplus,mminus);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppip);
      normalize(_h_K0pip);
      normalize(_h_K0Kp );
      normalize(_dalitz[0]);
      normalize(_h_Kppi0);
      normalize(_h_Kmpi0);
      normalize(_h_KpKm );
      normalize(_dalitz[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppip,_h_K0pip,_h_K0Kp;
    Histo1DPtr _h_Kppi0,_h_Kmpi0,_h_KpKm;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2015_I1403544);

}
