// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K+ KS0 pi0
  class BESIII_2021_I1859124 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1859124);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DP(ufs);
      DP.addStable(PID::PI0);
      DP.addStable(PID::K0S);
      declare(DP,"DP");
      // histos
      book(_h_Kppi0,1,1,1);
      book(_h_K0pi0,1,1,2);
      book(_h_KpK0,1,1,3);
      book(_dalitz, "dalitz",50,0.3,2.,50,0.3,2.);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{ 111,1}, {310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 111,1}, {310,1}};
      DecayedParticles DP = apply<DecayedParticles>(event, "DP");
      // loop over particles
      for(unsigned int ix=0;ix<DP.decaying().size();++ix) {
	int sign = 1;
	if (DP.decaying()[ix].pid()>0 && DP.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (DP.decaying()[ix].pid()<0 && DP.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & K0  = DP.decayProducts()[ix].at(      310)[0];
	const Particle & pi0 = DP.decayProducts()[ix].at(      111)[0];
	const Particle & Kp  = DP.decayProducts()[ix].at( sign*321)[0];
	double mplus  = (Kp.momentum() +pi0.momentum()).mass2();
	double mneut  = (K0.momentum() +pi0.momentum()).mass2();
	double mKK    = (Kp.momentum() + K0.momentum()).mass2();
	_h_Kppi0->fill(sqrt(mplus));
	_h_K0pi0->fill(sqrt(mneut));
	_h_KpK0 ->fill(sqrt(mKK  ));
	_dalitz ->fill(mplus,mneut);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Kppi0);
      normalize(_h_K0pi0);
      normalize(_h_KpK0 );
      normalize(_dalitz );
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_Kppi0,_h_K0pi0,_h_KpK0;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1859124);

}
