// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"

#ifndef RIVET_ENABLE_HEPMC_3
#include "HepMC/HepMCDefs.h"
#endif

namespace Rivet {


  /// Analysis of the generated cross-section
  class MC_XS : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(MC_XS);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      /// @todo Convert to Scatter1D or Counter
      book(_h_XS, "XS", {-0.5, 0.5});
      book(_h_N, "N", 1, 0.0, 1.0);
      book(_h_pmXS, "pmXS", 2, -1.0, 1.0);
      book(_h_pmN, "pmN", 2, -1.0, 1.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const size_t numWeights = event.weights().size();
      const vector<pair<double,double>> xsecs = event.crossSections();
      for (size_t m = 0; m < numWeights; ++m) {
        #if defined RIVET_ENABLE_HEPMC_3 || defined HEPMC_HAS_CROSS_SECTION
        size_t idx = (xsecs.size() == numWeights)? m : 0;
        const double xs    = xsecs[idx].first;
        const double xserr = xsecs[idx].second;
        _h_XS.get()->_getPersistent(m)->point(0).setY(xs, xserr);
        # endif
        const double weight = event.weights()[m];
        _h_pmXS.get()->_getPersistent(m)->fill(0.5*(weight > 0 ? 1. : -1), abs(weight));
        _h_pmN.get()->_getPersistent(m)->fill(0.5*(weight > 0 ? 1. : -1), 1.);
        _h_N.get()->_getPersistent(m)->fill(0.5, 1.0);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pmXS, crossSection()/sumOfWeights());
      #if !defined RIVET_ENABLE_HEPMC_3 && !defined HEPMC_HAS_CROSS_SECTION
      _h_XS->point(0).setY(crossSection(), 0.);
      #endif
    }

    /// @}


    /// @name Histograms
    /// @{
    Scatter2DPtr _h_XS;
    Histo1DPtr _h_pmXS, _h_pmN, _h_N;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(MC_XS);

}
