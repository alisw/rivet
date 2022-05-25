// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D* +/- production at 29 GeV
  class MARKII_1982_I177606 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1982_I177606);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      //Histograms
      book(_h_spect[0],2,1,1);
      book(_h_spect[1],3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==413)) {
	const double xp = 2.*p.E()/sqrtS();
	const double beta = p.p3().mod() / p.E();
	_h_spect[0]->fill(xp,1./beta);
	_h_spect[1]->fill(xp,1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spect[0], sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_h_spect[1], sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_spect[2];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1982_I177606);

}
