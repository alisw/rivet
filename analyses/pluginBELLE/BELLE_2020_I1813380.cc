// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c+ -> eta Lambda pi+
  class BELLE_2020_I1813380 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2020_I1813380);


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
      LAMBDAC.addStable(PID::LAMBDA);
      LAMBDAC.addStable(-PID::LAMBDA);
      declare(LAMBDAC, "LAMBDAC");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { PID::LAMBDA,1}, { 211,1}, { 221,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-PID::LAMBDA,1}, {-211,1}, { 221,1}};
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
	const Particle & eta = LAMBDAC.decayProducts()[ix].at( 221)[0];
	const Particle & pip  = LAMBDAC.decayProducts()[ix].at( sign*211)[0];
	_h[0]->fill((lam.momentum()+eta.momentum()).mass());
	_h[1]->fill((lam.momentum()+pip.momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BELLE_2020_I1813380);

}
