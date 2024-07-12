// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D*0 production at 34.4 GeV
  class JADE_1984_I221004 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1984_I221004);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_x    , 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==423)) {
	double xE =2.*p.E()/sqrtS();
	_h_x->fill(xE);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_x, crossSection()/microbarn/sumOfWeights()*sqr(sqrtS()));
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(JADE_1984_I221004);

}
