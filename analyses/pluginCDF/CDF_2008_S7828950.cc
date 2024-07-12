// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF Run II inclusive jet cross-section using the Midpoint algorithm.
  ///
  /// The analysis includes 1.1fb^-1 of CDF data and is the first with a
  /// cone algorithm to include the forward region of the detector.
  /// arXiv:0807.2204 to be published in PRD
  class CDF_2008_S7828950 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2008_S7828950);


    /// @name Analysis methods
    //@{

    // Book histos and set counters for number of events passed in each one
    void init() {
      const FinalState fs;
      declare(FastJets(fs, FastJets::CDFMIDPOINT, 0.7), "JetsM07");

      {Histo1DPtr tmp; _binnedHistosR07.add(  0, 0.1, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosR07.add(0.1, 0.7, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosR07.add(0.7, 1.1, book(tmp, 3, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosR07.add(1.1, 1.6, book(tmp, 4, 1, 1));}
      {Histo1DPtr tmp; _binnedHistosR07.add(1.6, 2.1, book(tmp, 5, 1, 1));}
    }


    // Do the analysis
    void analyze(const Event& event) {
      for (const Jet& jet : apply<FastJets>(event, "JetsM07").jets(Cuts::pT > 62*GeV)) {
        _binnedHistosR07.fill(jet.absrap(), jet.pT(), 1.0);
      }
    }


    // Normalise histograms to cross-section
    void finalize() {
      _binnedHistosR07.scale(crossSection()/nanobarn/sumOfWeights()/2.0, this);
    }

    //@}


  private:

    /// Histograms in different eta regions
    BinnedHistogram _binnedHistosR07;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2008_S7828950, CDF_2008_I790693);

}
