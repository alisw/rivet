// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> Lambda Sigmabar+- pi-+
  class BESIII_2013_I1261765 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2013_I1261765);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable( PID::PI0);
      psi.addStable( PID::K0S);
      psi.addStable( PID::SIGMAPLUS);
      psi.addStable( PID::SIGMAMINUS);
      psi.addStable(-PID::SIGMAPLUS);
      psi.addStable(-PID::SIGMAMINUS);
      psi.addStable( PID::LAMBDA);
      psi.addStable(-PID::LAMBDA);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 3122,1}, {-3222,1}, { 211,1} };
      static const map<PdgId,unsigned int> & mode1CC = { {-3122,1}, { 3222,1}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2   = { { 3122,1}, {-3112,1}, {-211,1} };
      static const map<PdgId,unsigned int> & mode2CC = { {-3122,1}, { 3112,1}, { 211,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	int ipi=211, isigma=-3222,sign=1;
	if(psi.modeMatches(ix,3,mode1)) {
	  sign =  1;
	}
	else if(psi.modeMatches(ix,3,mode1CC)) {
	  sign = -1;
	}
	else if(psi.modeMatches(ix,3,mode2)) {
	  ipi=-211;
	  isigma=-3112;
	  sign =  1;
	}
	else if(psi.modeMatches(ix,3,mode2CC)) {
	  ipi=-211;
	  isigma=-3112;
	  sign = -1;
	} 
	else
	  continue;
	const Particle & lam   = psi.decayProducts()[ix].at( sign*3122  )[0];
	const Particle & pim   = psi.decayProducts()[ix].at( sign*ipi   )[0];
	const Particle & sigma = psi.decayProducts()[ix].at( sign*isigma)[0];
	_h[0]->fill((lam.momentum()+pim  .momentum()).mass());
	_h[1]->fill((pim.momentum()+sigma.momentum()).mass());
	_h[2]->fill((lam.momentum()+sigma.momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BESIII_2013_I1261765);

}
