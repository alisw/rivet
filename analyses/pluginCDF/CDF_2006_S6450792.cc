// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF inclusive-jet cross-section differential in \f$ p_\perp \f$
  class CDF_2006_S6450792 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2006_S6450792);


    /// @name Analysis methods
    //@{

    void init() {
      FinalState fs;
      declare(FastJets(fs, FastJets::CDFMIDPOINT, 0.7), "ConeFinder");
      book(_h_jet_pt ,1, 1, 1);
    }


    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ConeFinder").jets(Cuts::pT > 61*GeV);
      for (const Jet& jet : jets) {
        if (inRange(jet.absrap(), 0.1, 0.7))
          _h_jet_pt->fill(jet.pT()/GeV);
      }
    }


    void finalize() {
      const double delta_y = 1.2;
      scale(_h_jet_pt, crossSection()/nanobarn/sumOfWeights()/delta_y);
    }

    //@}


  private:

    /// Histogram
    Histo1DPtr _h_jet_pt;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2006_S6450792, CDF_2006_I699933);

}
