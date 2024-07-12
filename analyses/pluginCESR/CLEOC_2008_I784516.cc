// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  psi(2S) -> J/Psi pi pi
  class CLEOC_2008_I784516 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2008_I784516);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::JPSI);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1 = { { 211,1}, {-211,1}, { 443,1} };
      static const map<PdgId,unsigned int> & mode2 = { { 111,2}          , { 443,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(psi.modeMatches(ix,3,mode1)) {
	  const Particle & pip = psi.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = psi.decayProducts()[ix].at(-211)[0];
	  _h[0]->fill((pip.momentum()+pim.momentum()).mass());
	}
	else if(psi.modeMatches(ix,3,mode2)) {
	  const Particles & pi0 = psi.decayProducts()[ix].at( 111);
	  _h[1]->fill((pi0[0].momentum()+pi0[1].momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)  normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2008_I784516);

}
