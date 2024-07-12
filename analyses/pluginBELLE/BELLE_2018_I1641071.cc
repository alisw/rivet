// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Omega_c decays
  class BELLE_2018_I1641071 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1641071);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4332);
      declare(ufs, "UFS");
      DecayedParticles OMEGAC(ufs);
      OMEGAC.addStable( PID::PI0);
      OMEGAC.addStable( PID::K0S);
      OMEGAC.addStable( PID::XIMINUS);
      OMEGAC.addStable(-PID::XIMINUS);
      OMEGAC.addStable( PID::XI0);
      OMEGAC.addStable(-PID::XI0);
      OMEGAC.addStable( PID::OMEGAMINUS);
      OMEGAC.addStable(-PID::OMEGAMINUS);
      declare(OMEGAC, "OMEGAC");
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 3334,1}, { 211,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode1CC = { {-3334,1}, {-211,1}, { 111,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 3312,1}, {-321,1}, { 211,2} };
      static const map<PdgId,unsigned int> & mode2CC = { {-3312,1}, { 321,1}, {-211,2} };
      static const map<PdgId,unsigned int> & mode3   = { { 3322,1}, {-321,1}, { 211,1} };
      static const map<PdgId,unsigned int> & mode3CC = { {-3322,1}, { 321,1}, {-211,1} };
      DecayedParticles OMEGAC = apply<DecayedParticles>(event, "OMEGAC");
      // loop over particles
      for(unsigned int ix=0;ix<OMEGAC.decaying().size();++ix) {
	int sign = OMEGAC.decaying()[ix].pid()/OMEGAC.decaying()[ix].abspid();
	// omega pi+ pi0
	if (( sign==1  && OMEGAC.modeMatches(ix,3,mode1  )) ||
	    ( sign==-1 && OMEGAC.modeMatches(ix,3,mode1CC))) {
	  const Particle & pi0   = OMEGAC.decayProducts()[ix].at( 111)[0];
	  const Particle & pip   = OMEGAC.decayProducts()[ix].at( sign*211 )[0];
	  _h[0]->fill((pi0.momentum()+pip.momentum()).mass());
	}
	// Xi- K- pi+pi+
	else if (( sign==1  && OMEGAC.modeMatches(ix,4,mode2  )) ||
		 ( sign==-1 && OMEGAC.modeMatches(ix,4,mode2CC))) {
	  const Particles & pip = OMEGAC.decayProducts()[ix].at( sign*211);
	  const Particle  & Km  = OMEGAC.decayProducts()[ix].at(-sign*321)[0];
	  const Particle  & xim = OMEGAC.decayProducts()[ix].at(sign*3312)[0];
	  _h[1]->fill((xim.momentum()+pip[0].momentum()).mass());
	  _h[1]->fill((xim.momentum()+pip[1].momentum()).mass());
	  _h[2]->fill(( Km.momentum()+pip[0].momentum()).mass());
	  _h[2]->fill(( Km.momentum()+pip[1].momentum()).mass());
	}
	// Xi0 K- pi+
	else if (( sign==1  && OMEGAC.modeMatches(ix,3,mode3  )) ||
		 ( sign==-1 && OMEGAC.modeMatches(ix,3,mode3CC))) {
	  const Particle & pip = OMEGAC.decayProducts()[ix].at( sign*211)[0];
	  const Particle & Km  = OMEGAC.decayProducts()[ix].at(-sign*321)[0];
	  _h[3]->fill(( Km.momentum()+pip.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1641071);

}
