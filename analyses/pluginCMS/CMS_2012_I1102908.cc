// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <sstream>

namespace Rivet {


  /// @brief CMS inclusive and exclusive dijet production ratio at large rapidity intervals
  class CMS_2012_I1102908 : public Analysis {
  public:

    CMS_2012_I1102908()
      : Analysis("CMS_2012_I1102908")
    {  }


  void init() {
    // Projections
    declare(FastJets(FinalState(), FastJets::ANTIKT, 0.5), "antikT");

    // Histograms
    /// @todo Can we manage to only register these as they are "really" created in the finalize()?
    _h_dijet_ratio = bookScatter2D(1, 1, 1);
    _h_MN_dijet_ratio = bookScatter2D(2, 1, 1);

    // Temporary histograms (directly instantiated)
    _h_DeltaY_exclusive = bookHisto1D("TMP/excl",refData(1, 1, 1));
    _h_DeltaY_inclusive = bookHisto1D("TMP/incl",refData(1, 1, 1));
    _h_DeltaY_MN = bookHisto1D("TMP/YMN",refData(1, 1, 1));
  }


  void analyze(const Event & event) {
    const double weight = event.weight();

    // Jets with  pT > 35.0, -4.7 < y < 4.7
    const JetAlg& jet_alg = apply<JetAlg>(event, "antikT");
    const Jets& jets = jet_alg.jets(Cuts::pT > 35*GeV && Cuts::absrap < 4.7);

    // Veto event if number of jets less than 2
    if (jets.size() < 2) return;

    // Loop over jet pairs
    double deltaY_MN = 0.0;
    for (size_t ij1 = 0; ij1 < jets.size(); ++ij1) {
      for (size_t ij2 = ij1 + 1; ij2 < jets.size(); ++ij2) {
        const double deltaY = fabs(jets[ij1].rapidity() - jets[ij2].rapidity());
        // Exclusive dijet case:
        if (jets.size() == 2) _h_DeltaY_exclusive->fill(deltaY, weight);
        // Inclusive jets case:
        _h_DeltaY_inclusive->fill(deltaY, weight);
        // Mueller-Navelet:
        if (deltaY > deltaY_MN) deltaY_MN = deltaY;
      }
    }
    _h_DeltaY_MN->fill(deltaY_MN, weight);
  }


  void finalize() {
    efficiency(_h_DeltaY_exclusive, _h_DeltaY_inclusive, _h_dijet_ratio);
    efficiency(_h_DeltaY_exclusive, _h_DeltaY_MN, _h_MN_dijet_ratio);
    transformY(*_h_dijet_ratio, _invert);
    transformY(*_h_MN_dijet_ratio, _invert);
  }


  private:

    /// Reciprocal function with div-by-zero protection, for inverting the efficiency measure
    static double _invert(double x) { return (x > 0) ? 1/x : 0; }

    /// @name Histograms
    //@{
    Scatter2DPtr _h_dijet_ratio, _h_MN_dijet_ratio;
    Histo1DPtr _h_DeltaY_inclusive, _h_DeltaY_exclusive, _h_DeltaY_MN;
    //@}

  };


  DECLARE_RIVET_PLUGIN(CMS_2012_I1102908);

}
