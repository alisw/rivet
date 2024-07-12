// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief D* production at 36.2 GeV
  class TASSO_1989_I278856 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1989_I278856);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_x_A1, 1, 1, 1);
      book(_h_x_A2, 1, 1, 2);
      book(_h_x_B1, 2, 1, 1);
      book(_h_x_B2, 2, 1, 2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==413)) {
	double modp = p.p3().mod();
	double xE = p.E()/meanBeamMom;
	double beta = modp/p.E();
	_h_x_A1->fill(xE);
	_h_x_A2->fill(xE,1./beta);
	_h_x_B1->fill(xE);
	_h_x_B2->fill(xE,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_x_A1, crossSection()/picobarn/sumOfWeights());
      scale(_h_x_A2, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale(_h_x_B1, crossSection()/picobarn/sumOfWeights());
      scale(_h_x_B2, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_x_A1, _h_x_A2, _h_x_B1, _h_x_B2;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TASSO_1989_I278856);


}
