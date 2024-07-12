// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DecayedParticles.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta' -> pi+pi- gamma decays
  class CRYSTAL_BARREL_1997_I456942 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BARREL_1997_I456942);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==331);
      declare(ufs, "UFS");
      DecayedParticles ETA(ufs);
      ETA.addStable(PID::PI0);
      declare(ETA, "ETA");
      // Book histograms
      book(_h_m, 1, 1, 1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { {211,1}, {-211,1}, { 22,1} };
      // Loop over eta' mesons
      DecayedParticles ETA = apply<DecayedParticles>(event, "ETA");
      // loop over particles
      for(unsigned int ix=0;ix<ETA.decaying().size();++ix) {
	// select right decay mode
	if ( !ETA.modeMatches(ix,3,mode)) continue;
	const Particle & pip = ETA.decayProducts()[ix].at( 211)[0];
	const Particle & pim = ETA.decayProducts()[ix].at(-211)[0];
	_h_m->fill((pip.momentum()+pim.momentum()).mass()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_m,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_m;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BARREL_1997_I456942);

}
