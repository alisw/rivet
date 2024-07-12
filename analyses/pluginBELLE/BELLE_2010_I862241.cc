// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> K+ K- KS0
  class BELLE_2010_I862241 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2010_I862241);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::K0S);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<4;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 321,1}, {-321,1}, { 310,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
	if (!B0.modeMatches(ix,3,mode)) continue;
	int sign = B0.decaying()[ix].pid()>0 ? 1 : -1;
      	const Particle & Kp  = B0.decayProducts()[ix].at( 321*sign)[0];
      	const Particle & Km  = B0.decayProducts()[ix].at(-321*sign)[0];
      	const Particle & K0  = B0.decayProducts()[ix].at( 310     )[0];
	_h[0]->fill((K0.momentum()+Kp.momentum()).mass());
	_h[1]->fill((K0.momentum()+Km.momentum()).mass());
	_h[2]->fill((Km.momentum()+Kp.momentum()).mass());
	_h[3]->fill((Km.momentum()+Kp.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix],1.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2010_I862241);

}
