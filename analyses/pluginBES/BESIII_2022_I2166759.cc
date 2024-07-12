// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/.psi anbd psi(2S) -> eta Sigma+ Sigmabar-
  class BESIII_2022_I2166759 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2166759);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443 or
						Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable( PID::ETA);
      psi.addStable( PID::SIGMAPLUS);
      psi.addStable( PID::SIGMAMINUS);
      psi.addStable(-PID::SIGMAPLUS);
      psi.addStable(-PID::SIGMAMINUS);
      declare(psi, "psi");
      // histos
      for(unsigned int ix=0;ix<6;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 3222,1}, {-3222,1}, { 221,1} };
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	unsigned int ioff = psi.decaying()[ix].pid()==443 ?  0 : 3;
	const Particle & sigp = psi.decayProducts()[ix].at( 3222)[0];
	const Particle & sigm = psi.decayProducts()[ix].at(-3222)[0];
	const Particle & eta  = psi.decayProducts()[ix].at(  221)[0];
	_h[ioff+0]->fill((sigp.momentum()+eta .momentum()).mass());
	_h[ioff+1]->fill((sigm.momentum()+eta .momentum()).mass());
	_h[ioff+2]->fill((sigp.momentum()+sigm.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<6;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2166759);

}
