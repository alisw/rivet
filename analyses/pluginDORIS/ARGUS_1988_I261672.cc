// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Sigma_c 0,++ spectrum
  class ARGUS_1988_I261672 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1988_I261672);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_x, 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4112 or Cuts::abspid==4222)) {
	const double xp = 2.*p.p3().mod()/sqrtS();
	_h_x->fill(xp);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ARGUS_1988_I261672);


}
