// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1982_I179022 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1982_I179022);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_rho, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==113)) {
	double xE = p.E()/meanBeamMom;
	double beta = p.p3().mod()/p.E();
	_h_rho->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_rho, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights()); // norm to cross section

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_rho;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1982_I179022);


}
