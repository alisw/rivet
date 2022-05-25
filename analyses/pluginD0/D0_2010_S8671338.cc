// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Measurement of Z(->muon muon) pT differential cross-section
  ///
  /// @author Flavia Dias
  class D0_2010_S8671338 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_S8671338);


    ///@name Analysis methods
    //@{

    /// Add projections and book histograms
    void init() {
      Cut cut = Cuts::abseta < 1.7 && Cuts::pT > 15*GeV;
      ZFinder zfinder(FinalState(), cut, PID::MUON, 65*GeV, 115*GeV, 0.2, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");

      book(_h_Z_pT_normalised ,1, 1, 1);
      book(_h_Z_pT_xs ,2, 1, 1);
    }


    // Do the analysis
    void analyze(const Event& e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size()==1) {
        double ZpT = zfinder.bosons()[0].pT()/GeV;
        _h_Z_pT_normalised->fill(ZpT);
        _h_Z_pT_xs->fill(ZpT);
      }
    }


    /// Finalize
    void finalize() {
      normalize(_h_Z_pT_normalised);
      scale(_h_Z_pT_xs, crossSection()/sumOfWeights());
    }

    //@}


  private:

    /// @name Histogram
    Histo1DPtr _h_Z_pT_normalised;
    Histo1DPtr _h_Z_pT_xs;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_S8671338, D0_2010_I856972);

}
