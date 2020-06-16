// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"

#ifndef RIVET_ENABLE_HEPMC_3
#include "HepMC/HepMCDefs.h"
#endif

namespace Rivet {

  /// @brief Analysis for the generated cross section
  class MC_XS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_XS);

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      /// @todo Convert to Scatter1D or Counter
      book(_h_XS, "XS");
      book(_h_N, "N", 1, 0.0, 1.0);
      book(_h_pmXS, "pmXS", 2, -1.0, 1.0);
      book(_h_pmN, "pmN", 2, -1.0, 1.0);

      _mc_xs = _mc_error = 0.;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      #if defined RIVET_ENABLE_HEPMC_3
      //@todo HepMC3::GenCrossSection methods aren't const accessible :(
      RivetHepMC::GenCrossSection gcs = *(event.genEvent()->cross_section());
      _mc_xs    = gcs.xsec();
      _mc_error = gcs.xsec_err();
      #elif defined HEPMC_HAS_CROSS_SECTION
      _mc_xs    = event.genEvent()->cross_section()->cross_section();
      _mc_error = event.genEvent()->cross_section()->cross_section_error();
      #endif // VERSION_CODE >= 3000000

      const size_t numWeights = handler().numWeights();
      const vector<size_t>& indices = handler().weightIndices();
      assert(numWeights == indices.size());
      for (size_t m = 0; m < numWeights; ++m) {
        const double weight = event.weights()[indices[m]];
        _h_pmXS.get()->_getPersistent(m)->fill(0.5*(weight > 0 ? 1. : -1), abs(weight));
        _h_pmN.get()->_getPersistent(m)->fill(0.5*(weight > 0 ? 1. : -1), 1.);
        _h_N.get()->_getPersistent(m)->fill(0.5, 1.0);
      }
      
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pmXS, crossSection()/sumOfWeights());
      #ifndef HEPMC_HAS_CROSS_SECTION
      _mc_xs = crossSection();
      _mc_error = 0.0;
      #endif
      _h_XS->addPoint(0, _mc_xs, 0.5, _mc_error);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Scatter2DPtr _h_XS;
    Histo1DPtr _h_pmXS, _h_pmN, _h_N;
    double _mc_xs, _mc_error;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_XS);

}
