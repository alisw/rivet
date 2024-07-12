// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> p 2K0S and p K0S eta
  class BELLE_2022_I2160445 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2160445);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::K0S);
      LAMBDAC.addStable(PID::ETA);
      declare(LAMBDAC, "LAMBDAC");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2+ix;++iy)
	  book(_h[ix][iy],1+ix,1,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { PID::PROTON,1}, {310,2}};
      static const map<PdgId,unsigned int> & mode1CC = { {-PID::PROTON,1}, {310,2}};
      static const map<PdgId,unsigned int> & mode2   = { { PID::PROTON,1}, {310,1}, { 221,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-PID::PROTON,1}, {310,1}, { 221,1}};
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for(unsigned int ix=0;ix<LAMBDAC.decaying().size();++ix) {
	int sign = 1, mode=-1;
	if (LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode1)) {
	  sign=1;
	  mode=0;
	}
	else if  (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,mode1CC)) {
	  sign=-1;
	  mode=0;
	}
	else if (LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode2)) {
	  sign=1;
	  mode=1;
	}
	else if  (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,mode2CC)) {
	  sign=-1;
	  mode=1;
	} 
	else
	  continue;
	const Particle  & pp  = LAMBDAC.decayProducts()[ix].at( sign*PID::PROTON)[0];
	const Particles & K0  = LAMBDAC.decayProducts()[ix].at( 310);
	if(mode==0) {
	  _h[0][0]->fill((pp   .momentum()+K0[0].momentum()).mass2());
	  _h[0][0]->fill((pp   .momentum()+K0[1].momentum()).mass2());
	  _h[0][1]->fill((K0[0].momentum()+K0[1].momentum()).mass2());
	}
	else {
	  const Particle & eta = LAMBDAC.decayProducts()[ix].at( 221)[0];
	  _h[1][0]->fill((K0[0].momentum()+eta.momentum()).mass2());
	  _h[1][1]->fill((pp.momentum()+ K0[0].momentum()).mass2());
	  _h[1][2]->fill((pp.momentum()+eta.momentum()).mass2());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2+ix;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2160445);

}
