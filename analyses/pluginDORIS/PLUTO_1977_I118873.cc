// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Kaon spectra at 3.63, 4.03 and 4.5 GeV
  class PLUTO_1977_I118873 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PLUTO_1977_I118873);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      
      if      (fuzzyEquals(sqrtS()/GeV, 3.63, 1E-3)) {
	book(_h_spectrum, 2, 1, 1);
      }
      else if (fuzzyEquals(sqrtS()/GeV, 4.03, 1E-3)) {
	book(_h_spectrum, 3, 1, 1);
      }
      else if (fuzzyEquals(sqrtS()/GeV, 4.5, 1E-3)) {
	book(_h_spectrum, 4, 1, 1);
      }
      else
	MSG_ERROR("Beam energy not supported!");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // unstable particles
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").
	       particles(Cuts::pid==PID::K0S)) {
	double xp = p.E()/meanBeamMom;
	double beta = p.p3().mod()/p.E();
	_h_spectrum->fill(xp,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spectrum, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_spectrum;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PLUTO_1977_I118873);


}
