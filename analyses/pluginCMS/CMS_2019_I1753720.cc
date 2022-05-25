// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ttbb cross section all-jet 2016 data
  class CMS_2019_I1753720 : public Analysis {
    public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_I1753720);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Jets
      // Only use visible particles with |eta|<6 (miniAOD content)
      declare(FastJets(VisibleFinalState(Cuts::abseta < 6.), FastJets::ANTIKT, 0.4), "Jets");

      // Book xsec histo
      book(_hist_xsec_fid, "d01-x01-y01");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
       const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 20*GeV);
       const Jets jets_30 = filter_select(jets, [](const Jet& j) { return j.pT() > 30*GeV; } );
       const Jets bjets = filter_select(jets, [](const Jet& j) { return j.bTagged(); } );

       if (jets.size() >= 8 && jets_30.size() >= 6 && bjets.size() >= 4) {
           _hist_xsec_fid->fill(1.);
       }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_hist_xsec_fid, crossSection()/picobarn/sumOfWeights());

    }

    //@}
    
    private:

    /// Histogram for fiducial cross section
    Histo1DPtr _hist_xsec_fid;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2019_I1753720);


}
