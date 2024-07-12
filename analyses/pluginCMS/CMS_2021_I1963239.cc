// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief  Measurement of inclusive and Mueller-Navelet dijet cross sections and their ratios at 2.76 TeV
  class CMS_2021_I1963239 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1963239);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 5.2);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.5
      FastJets jetfs(fs, FastJets::ANTIKT, 0.5);
      declare(jetfs, "jets");

      // Book histograms
      // specify custom binning
      // take binning from reference data using HEPData ID (digits in "d01-x01-y01" etc.)
      
      book(_h["inclusive"], 7, 1, 1);
      book(_h["MN"], 8, 1, 1);
      book(_s["R_incl"], 9, 1, 1);
      book(_s["R_incl_veto"], 10, 1, 1);
      book(_s["R_MN"], 11, 1, 1);
      book(_s["R_MN_veto"], 12, 1, 1);

      // Temporary histograms (directly instantiated)
      book(_h["exclusive"], "_exclusive", refData(7, 1, 1));
      book(_h["exclusive_veto"], "_exclusive_veto", refData(7, 1, 1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets20 = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::absrap < 4.7);
      Jets jets35 = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 35*GeV && Cuts::absrap < 4.7);

      if (jets35.size() < 2) return;

      // Loop over jet pairs
      double deltaY_MN = 0.0;
      for (size_t ij1 = 0; ij1 < jets35.size(); ++ij1) {
        for (size_t ij2 = ij1 + 1; ij2 < jets35.size(); ++ij2) {
          const double deltaY = fabs(jets35[ij1].rapidity() - jets35[ij2].rapidity());
          // Exclusive dijet case:
          if (jets35.size() == 2) {
            _h["exclusive"]->fill(deltaY, weight);
            //Exclusive with veto 20 GeV dijet case:
            if (jets20.size() == 2) {
              _h["exclusive_veto"]->fill(deltaY, weight);
            }
          }
          // Inclusive jets case:
          _h["inclusive"]->fill(deltaY, weight);
          // Mueller-Navelet:
          if (deltaY > deltaY_MN) deltaY_MN = deltaY;
        }
      }
      // Fill histogram with MN dijets Delta y
      _h["MN"]->fill(deltaY_MN, weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Calculate ratios
      efficiency(_h["exclusive"], _h["inclusive"], _s["R_incl"]);
      efficiency(_h["exclusive"], _h["MN"], _s["R_MN"]);
      efficiency(_h["exclusive_veto"], _h["inclusive"], _s["R_incl_veto"]);
      efficiency(_h["exclusive_veto"], _h["MN"], _s["R_MN_veto"]);

      transformY(*_s["R_incl"], _invert);
      transformY(*_s["R_MN"], _invert);
      transformY(*_s["R_incl_veto"], _invert);
      transformY(*_s["R_MN_veto"], _invert);


      scale(_h["inclusive"], crossSection()/picobarn/sumOfWeights()); // norm to generated cross-section in pb
      scale(_h["MN"], crossSection()/picobarn/sumOfWeights());        // norm to generated cross-section in pb

    }

    ///@}


    /// @name Histograms
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr>_s;
    ///@}
    private:
    
    /// Reciprocal function with div-by-zero protection, for inverting the efficiency measure
    static double _invert(double x) { return (x > 0) ? 1/x : 0; }


  };


  RIVET_DECLARE_PLUGIN(CMS_2021_I1963239);

}
