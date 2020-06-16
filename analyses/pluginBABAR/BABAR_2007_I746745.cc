// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Spectrum for Omega_c
  class BABAR_2007_I746745 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2007_I746745);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(),"UFS");
      book(_h_p,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idOmega = 4332;
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idOmega)) {
	_h_p->fill(p.momentum().p3().mod());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_p);
    }
    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_p;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2007_I746745);

}
