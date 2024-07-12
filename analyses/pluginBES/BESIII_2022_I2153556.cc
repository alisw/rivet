// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> Lambda pi+ pi0
  class BESIII_2022_I2153556 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2153556);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::LAMBDA);
      LAMBDAC.addStable(-PID::LAMBDA);
      declare(LAMBDAC, "LAMBDAC");
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { PID::LAMBDA,1}, { 211,1}, { 111,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-PID::LAMBDA,1}, {-211,1}, { 111,1}};
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for(unsigned int ix=0;ix<LAMBDAC.decaying().size();++ix) {
	int sign = 1;
	if (LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode)) {
	  sign=1;
	}
	else if  (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,modeCC)) {
	  sign=-1;
	}
	else
	  continue;
	const Particle & lam = LAMBDAC.decayProducts()[ix].at( sign*PID::LAMBDA)[0];
	const Particle & pi0 = LAMBDAC.decayProducts()[ix].at(      PID::PI0   )[0];
	const Particle & pip = LAMBDAC.decayProducts()[ix].at( sign*PID::PIPLUS)[0];
	_h[0]->fill((pip.momentum()+pi0.momentum()).mass());
	_h[1]->fill((lam.momentum()+pip.momentum()).mass());
	_h[2]->fill((lam.momentum()+pi0.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2153556);

}
