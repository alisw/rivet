// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief eta' -> gamma e+e-
  class BESIII_2015_I1364494 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1364494);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==331);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      declare(ETA, "ETA");
      // Book histograms
      book(_h_m, 1, 1, 3);
      book(_netap, "TMP/netap");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode0  = { { 22,2} };
      static const map<PdgId,unsigned int> & mode1  = { {11,1}, {-11,1}, { 22,1} };
      // Loop over eta' mesons
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	// refewrnece mode for denominator
	if(ETA.modeMatches(ix,2,mode0))
	   _netap->fill();
	// select right decay mode
	else if ( ETA.modeMatches(ix,3,mode1)) {
	  const Particle & ep = ETA.decayProducts()[ix].at( 11)[0];
	  const Particle & em = ETA.decayProducts()[ix].at(-11)[0];
	  _h_m->fill( (em.momentum()+ep.momentum()).mass());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // divide by no so BR and mult by bin width
      // and 100 as in %
      scale(_h_m,0.1*100./ *_netap);

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_m;
    CounterPtr _netap;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2015_I1364494);


}
