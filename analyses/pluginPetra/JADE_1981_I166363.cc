// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pbar and lambdabar at 34 GeV
  class JADE_1981_I166363 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1981_I166363);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_pbar     , 1, 1, 1);
      book(_h_lambdabar, 2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // at least 5 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      if (numParticles < 3) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==-2212 or Cuts::pid==-3122)) {
	if(p.pid()==-2212)
	  _h_pbar->fill(p.p3().mod());
	else
	  _h_lambdabar->fill(p.p3().mod());
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pbar     , crossSection()/nanobarn/sumOfWeights());
      scale(_h_lambdabar, crossSection()/nanobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr  _h_pbar, _h_lambdabar;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(JADE_1981_I166363);


}
