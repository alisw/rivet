// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Xi_c+ -> Sigma+ K- pi+ and Xi- pi+ pi+
  class FOCUS_2003_I618864 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FOCUS_2003_I618864);


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
      XICP.addStable(3312);
      XICP.addStable(-3312);
      declare(XICP, "XICP");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1+ix,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { PID::KMINUS,1}, { 3222,1}, { PID::PIPLUS ,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { PID::KPLUS ,1}, {-3222,1}, { PID::PIMINUS,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 3312,1}, { PID::PIPLUS ,2}};
      static const map<PdgId,unsigned int> & mode2CC = { {-3312,1}, { PID::PIMINUS,2}};
      DecayedParticles XICP = apply<DecayedParticles>(event, "XICP");
      // loop over particles
      for(unsigned int ix=0;ix<XICP.decaying().size();++ix) {
	int sign = XICP.decaying()[ix].pid()/XICP.decaying()[ix].abspid();
	if ( (sign== 1 && XICP.modeMatches(ix,3,mode1  )) ||
	     (sign==-1 && XICP.modeMatches(ix,3,mode1CC)) ) {
	  const Particle & pip = XICP.decayProducts()[ix].at( sign*PID::PIPLUS)[0];
	  const Particle & Km  = XICP.decayProducts()[ix].at( sign*PID::KMINUS)[0];
	  for(unsigned int ix=0;ix<2;++ix)
	    _h[0]->fill((pip.momentum()+Km.momentum()).mass());
	}
	else if ( (sign== 1 && XICP.modeMatches(ix,3,mode2  )) ||
		  (sign==-1 && XICP.modeMatches(ix,3,mode2CC)) ) {
	  const Particles & pip = XICP.decayProducts()[ix].at( sign*PID::PIPLUS);
	  const Particle  & xi  = XICP.decayProducts()[ix].at( sign*3312)[0];
	  _h[1]->fill((pip[0].momentum()+xi.momentum()).mass());
	  _h[1]->fill((pip[1].momentum()+xi.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FOCUS_2003_I618864);

}
