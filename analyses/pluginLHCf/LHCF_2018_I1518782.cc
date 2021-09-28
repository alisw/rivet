// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {

  /// @brief Forward photon production cross-section at 13 TeV
  class LHCF_2018_I1518782 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(LHCF_2018_I1518782);

    /// @name Analysis methods
    //@{

    void init() {
      declare(FinalState(), "FS");
      book(_h_n_en_eta1, 1, 1, 1);
      book(_h_n_en_eta2, 2, 1, 1);
    }

    void analyze(const Event& event) {

      // Select photons above threshold
      Particles fs_particles = apply<FinalState> (event, "FS").particles(Cuts::abspid==PID::PHOTON && Cuts::E>=200/GeV && Cuts::abseta>8.81);

      for (const Particle& p : fs_particles) {
        // Double analysis efficiency with a two-sided LHCf --- NOTE: taken care of in finalize division by 2
        const double eta = abs(p.eta());
        const double energy = p.E()/GeV;
        if ( eta > 10.94 ) _h_n_en_eta1->fill(energy);
        if ( eta <  8.99 ) _h_n_en_eta2->fill(energy);
      }
    }

    void finalize() {
      // Norm to cross-section
      scale(_h_n_en_eta1, crossSection()/millibarn/sumOfWeights()/2.);
      scale(_h_n_en_eta2, crossSection()/millibarn/sumOfWeights()/2.);
    }
    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_n_en_eta1, _h_n_en_eta2;
    //@}
  };

  DECLARE_RIVET_PLUGIN(LHCF_2018_I1518782);
}
