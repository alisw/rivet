// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta' -> pi+pi- gamma decay
  class BESIII_2018_I1641075 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1641075);


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
      book(_h_m, 1, 1, 5);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { PID::PIPLUS,1}, { PID::PIMINUS ,1}, { PID::GAMMA,1}};
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	if (!ETA.modeMatches(ix,3,mode)) continue;
	const Particle & pip = ETA.decayProducts()[ix].at( PID::PIPLUS )[0];
	const Particle & pim = ETA.decayProducts()[ix].at( PID::PIMINUS)[0];
	_h_m->fill((pip.momentum()+pim.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_m);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_m;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BESIII_2018_I1641075);


}
