// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta' -> pions
  class BESIII_2017_I1469067 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1469067);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==PID::ETAPRIME);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      ETA.addStable(PID::K0S);
      ETA.addStable(PID::ETA);
      declare(ETA, "ETA");
      // histograms
      for(unsigned int ix=0;ix<4;++ix) {
	book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1  = { {211,1}, {-211,1}, {111,1} };
      static const map<PdgId,unsigned int> & mode2   = { {111,3} };
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	// select right decay mode
	if ( ETA.modeMatches(ix,3,mode1)) {
	  const Particle & pi0 = ETA.decayProducts()[ix].at( 111)[0];
	  const Particle & pip = ETA.decayProducts()[ix].at( 211)[0];
	  const Particle & pim = ETA.decayProducts()[ix].at(-211)[0];
	  _h[0]->fill((pip.momentum()+pim.momentum()).mass());
	  _h[1]->fill((pip.momentum()+pi0.momentum()).mass());
	  _h[2]->fill((pim.momentum()+pi0.momentum()).mass());
	}
	else if  ( ETA.modeMatches(ix,3,mode2)) {
	  const Particles & pi0 = ETA.decayProducts()[ix].at(111);
	  for(unsigned int ix=0;ix<3;++ix) {
	    for(unsigned int iy=ix+1;iy<3;++iy)
	      _h[3]->fill((pi0[ix].momentum()+pi0[iy].momentum()).mass());
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize the histograms
      for(unsigned int ix=0;ix<4;++ix)
	normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1469067);

}
