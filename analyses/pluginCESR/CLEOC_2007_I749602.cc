// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  D+ -> pi+pi+pi-
  class CLEOC_2007_I749602 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2007_I749602);


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
      DP.addStable(PID::ETA);
      DP.addStable(PID::ETAPRIME);
      declare(DP, "DP");
      // histos
      book(_h_pippim,1,1,1);
      book(_h_pippip,1,1,2);
      book(_dalitz, "dalitz",50,0.,1.8,50,0.0,3.1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 211,2},{-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 211,1},{-211,2}};
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
	const Particles & pip= DP.decayProducts()[ix].at( sign*211);
	const Particles & pim= DP.decayProducts()[ix].at(-sign*211);
	double m1 = (pim[0].momentum()+pip[0].momentum()).mass2();
	double m2 = (pim[0].momentum()+pip[1].momentum()).mass2();
	double m3 = (pip[0].momentum()+pip[1].momentum()).mass2();
	if(m1>m2) swap(m1,m2);
	_dalitz->fill(m1,m2);
	// K_S0 veto
	if(m1<0.2 || m1>0.3) {
	  _h_pippim->fill(m1);
	  _h_pippim->fill(m2);
	  _h_pippip->fill(m3);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_pippim);
      normalize(_h_pippip);
      normalize(_dalitz);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_pippim,_h_pippip;
    Histo2DPtr _dalitz;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2007_I749602);

}
