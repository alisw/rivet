// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Xi_c+ -> Sigma+ K- pi+
  class CLEOII_1996_I397787 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1996_I397787);


    /// @name Analysis methods
    /// @{
    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4232);
      declare(ufs, "UFS");
      DecayedParticles XICP(ufs);
      XICP.addStable(PID::PI0);
      XICP.addStable(PID::K0S);
      XICP.addStable(PID::ETA);
      XICP.addStable(3222);
      XICP.addStable(-3222);
      declare(XICP, "XICP");
      // histograms
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { PID::KMINUS,1}, { 3222,1}, { PID::PIPLUS ,1}};
      static const map<PdgId,unsigned int> & modeCC = { { PID::KPLUS ,1}, {-3222,1}, { PID::PIMINUS,1}};
      DecayedParticles XICP = apply<DecayedParticles>(event, "XICP");
      // loop over particles
      for(unsigned int ix=0;ix<XICP.decaying().size();++ix) {
	int sign = 1;
	if (XICP.decaying()[ix].pid()>0 && XICP.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (XICP.decaying()[ix].pid()<0 && XICP.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & pip = XICP.decayProducts()[ix].at( sign*PID::PIPLUS)[0];
	const Particle & Km  = XICP.decayProducts()[ix].at( sign*PID::KMINUS)[0];
	_h->fill((pip.momentum()+Km.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1996_I397787);

}
