// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D0 -> KS0 pi+pi- and K+K-
  class BABAR_2010_I853279 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2010_I853279);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==421);
      declare(ufs, "UFS");
      DecayedParticles D0(ufs);
      D0.addStable(PID::PI0);
      D0.addStable(PID::K0S);
      D0.addStable(PID::ETA);
      D0.addStable(PID::ETAPRIME);
      declare(D0, "D0");
      // Histograms
      book(_h_Kpim,1,1,1);
      book(_h_Kpip,1,1,2);
      book(_h_pipi,1,1,3);
      book(_dalitz[0], "dalitz1",50,0.3,3.2,50,0.3,3.2);
      book(_h_K0Km,1,1,4);
      book(_h_K0Kp,1,1,5);
      book(_h_KpKm,1,1,6);
      book(_dalitz[1], "dalitz2",50,0.9,1.9,50,0.9,1.9);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & modePi  = { { 310,1}, { 211,1},{-211,1}};
      static const map<PdgId,unsigned int> & modeK   = { { 310,1}, { 321,1},{-321,1}};
      DecayedParticles D0 = apply<DecayedParticles>(event, "D0");
      // loop over particles
      for(unsigned int ix=0;ix<D0.decaying().size();++ix) {
	int sign = D0.decaying()[ix].pid()/421;
	// KS0 pi+pi-
	if (D0.modeMatches(ix,3,modePi) ) {
	  const Particles & pip= D0.decayProducts()[ix].at( sign*211);
	  const Particles & pim= D0.decayProducts()[ix].at(-sign*211);
	  const Particles & K0 = D0.decayProducts()[ix].at( 310);
	  double mminus = (pim[0].momentum()+K0[0].momentum() ).mass2();
	  double mplus  = (pip[0].momentum()+K0[0].momentum() ).mass2();
	  double mpipi  = (pip[0].momentum()+pim[0].momentum()).mass2();
	  _h_Kpip->fill(mplus);
	  _h_Kpim->fill(mminus);
	  _h_pipi->fill(mpipi);
	  _dalitz[0]->fill(mminus,mplus); 
	}
	else if (D0.modeMatches(ix,3,modeK) ) {
	  const Particles & Kp = D0.decayProducts()[ix].at( sign*321);
	  const Particles & Km = D0.decayProducts()[ix].at(-sign*321);
	  const Particles & K0 = D0.decayProducts()[ix].at( 310);
	  double mminus = (Km[0].momentum()+K0[0].momentum() ).mass2();
	  double mplus  = (Kp[0].momentum()+K0[0].momentum() ).mass2();
	  double mKK  = (Kp[0].momentum()+Km[0].momentum()).mass2();
	  _h_K0Kp->fill(mplus);
	  _h_K0Km->fill(mminus);
	  _h_KpKm->fill(mKK);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kpim);
      normalize(_h_Kpip);
      normalize(_h_pipi);
      normalize(_dalitz[0]);
      normalize(_h_K0Km);
      normalize(_h_K0Kp);
      normalize(_h_KpKm);
      normalize(_dalitz[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kpim,_h_pipi,_h_Kpip;
    Histo1DPtr _h_K0Km,_h_KpKm,_h_K0Kp;
    Histo2DPtr _dalitz[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2010_I853279);

}
