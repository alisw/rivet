// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF inclusive jet cross-section using the \f$ k_\perp \f$ algorithm
  class CDF_2007_S7057202 : public Analysis {
  public:

    CDF_2007_S7057202()
      : Analysis("CDF_2007_S7057202")
    {    }


    void init() {
      // Set up projections
      const FinalState fs;
      declare(FastJets(fs, FastJets::KT, 0.5), "JetsD05");
      declare(FastJets(fs, FastJets::KT, 0.7), "JetsD07");
      declare(FastJets(fs, FastJets::KT, 1.0), "JetsD10");

      // Book histos
      {Histo1DPtr tmp; _binnedHistosD07.add(  0, 0.1, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosD07.add(0.1, 0.7, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosD07.add(0.7, 1.1, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosD07.add(1.1, 1.6, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosD07.add(1.6, 2.1, book(tmp, 5, 1, 1));}
      book(_histoD05 ,6, 1, 1);
      book(_histoD10 ,7, 1, 1);
    }


    void analyze(const Event& event) {
      const double weight = 1.0;

      for (const Jet& jet : apply<JetAlg>(event, "JetsD07").jets(Cuts::pT > 54*GeV))
        _binnedHistosD07.fill(jet.absrap(), jet.pT(), weight);

      for (const Jet& jet : apply<JetAlg>(event, "JetsD05").jets(Cuts::pT > 54*GeV))
        if (inRange(jet.absrap(), 0.1, 0.7)) _histoD05->fill(jet.pT(), weight);

      for (const Jet& jet : apply<JetAlg>(event, "JetsD10").jets(Cuts::pT > 54*GeV))
        if (inRange(jet.absrap(), 0.1, 0.7)) _histoD10->fill(jet.pT(), weight);
    }


    // Normalise histograms to cross-section
    void finalize() {
      const double xSec = crossSectionPerEvent()/nanobarn;

      scale(_histoD05, xSec);
      scale(_histoD10, xSec);
      // scale to xSec/yBinWidth and take into account the double yBinWidth due
      // to the absolute value of y
      _binnedHistosD07.scale(xSec/2.0, this);
    }

  private:

    BinnedHistogram _binnedHistosD07;

    // Single histogram for the \f$R=0.5\f$ \f$k_\perp\f$ jets
    Histo1DPtr _histoD05;

    // Single histogram for the \f$R=1.0\f$ \f$k_\perp\f$ jets
    Histo1DPtr _histoD10;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2007_S7057202);

}
