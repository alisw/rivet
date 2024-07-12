// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> eta' e+e-
  class BESIII_2018_I1664315 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1664315);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==PID::PSI2S);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::ETAPRIME);
      declare(psi, "PSI");
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix], 1, 1, ix+1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // define the decay mode
      static const map<PdgId,unsigned int> & mode   = { {PID::ETAPRIME,1},{ 11,1}, { -11,1}};
      DecayedParticles psi = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for(unsigned int ix=0;ix<psi.decaying().size();++ix) {
	if(!psi.modeMatches(ix,3,mode)) continue;
	const Particle & em = psi.decayProducts()[ix].at( 11)[0];
	const Particle & ep = psi.decayProducts()[ix].at(-11)[0];
	double mee = (ep.momentum()+em.momentum()).mass();
	for(unsigned int ix=0;ix<4;++ix) _h[ix]->fill(mee);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1664315);

}
