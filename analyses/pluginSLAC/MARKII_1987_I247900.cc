// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Omega production at 29 GeV
  class MARKII_1987_I247900 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1987_I247900);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      //Histograms
      book(_h_spect,1,1,1);
      book(_h_sigma,2,1,1);
      book(_h_rate ,2,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==3334)) {
	const double xp = 2.*p.E()/sqrtS();
	const double beta = p.p3().mod() / p.E();
	_h_spect->fill(xp,1./beta);
	_h_sigma->fill(sqrtS());
	_h_rate->fill(sqrtS());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_spect, sqr(sqrtS())*crossSection()/nanobarn/sumOfWeights());
      scale(_h_sigma, crossSection()/picobarn/sumOfWeights());
      scale(_h_rate , 1./sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_spect,_h_sigma,_h_rate;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1987_I247900);

}
