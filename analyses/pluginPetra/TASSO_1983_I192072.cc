// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Xi- spectrum at 34 GeV
  class TASSO_1983_I192072 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1983_I192072);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_spectrum, 2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==3312)) {
	double modp = p.p3().mod();
	double beta = modp/p.E();
	double xE = p.E()/meanBeamMom;
	_h_spectrum->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spectrum, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_spectrum;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1983_I192072);


}
