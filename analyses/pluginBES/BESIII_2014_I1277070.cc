// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K_0S pi+ pi0
  class BESIII_2014_I1277070 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2014_I1277070);


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
      book(_h_pipi  ,1,1,1);
      book(_h_KS0pi0,1,1,2);
      book(_h_KS0pip,1,1,3);
      book(_dalitz, "dalitz",50,0.3,3.1,50,0.,2.0);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,1},{ 111,1}, {310,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-211,1},{ 111,1}, {310,1}};
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
	const Particle & pip = DP.decayProducts()[ix].at( sign*211)[0];
	double mpipi   = (pi0.momentum()+pip.momentum()).mass2();
	double mKS0pip = (K0 .momentum()+pip.momentum()).mass2();
	double mKS0pi0 = (K0 .momentum()+pi0.momentum()).mass2();
	_h_pipi  ->fill(mpipi);
	_h_KS0pip->fill(mKS0pip);
	_h_KS0pi0->fill(mKS0pi0);
	_dalitz->fill(mKS0pi0,mpipi);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pipi);
      normalize(_h_KS0pi0);
      normalize(_h_KS0pip);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pipi,_h_KS0pi0, _h_KS0pip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2014_I1277070);

}
