// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta production at 29 GeV
  class MARKII_1988_I261194 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1988_I261194);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      //Histograms
      book(_h_spect,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==221)) {
	const double xp = 2.*p.E()/sqrtS();
	const double beta = p.p3().mod() / p.E();
	_h_spect->fill(xp,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spect, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_spect;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1988_I261194);

}
