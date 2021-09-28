// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief ATLAS azimuthal decorrelation with jet veto analysis
  /// @author James Robinson <james.robinson@cern.ch>
  class ATLAS_2014_I1307243 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2014_I1307243);


    /// Book histograms and initialise projections before the run
    void init() {

      _dy_max = 8; _years = { 2010, 2011 };
      _Qnoughts = { 20., 30., 40., 50., 60., 70., 80., 90., 100. };
      /// Initialise and register projections here
      FastJets fastJets(FinalState(), FastJets::ANTIKT, 0.6, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(fastJets, "AntiKt6JetsWithInvisibles");


      /// Book histograms
      for (const string& cat : { "inclusive", "gap" }) {
        const size_t offset = (cat == "gap") ? 1 : 0;

        // Temporary inclusive and gap histograms
        book(_aux_dy[cat], "_" + cat + "_dy", refData(1, 1, 1));
        book(_aux_pTbar[cat], "_" + cat + "_pTbar", refData(2, 1, 1));

        book(_h_C2C1_dy[cat],    7 + 4 * offset, 1, 1, false);
        book(_h_C2C1_pTbar[cat], 8 + 4 * offset, 1, 1, false);

        // Azimuthal moment histograms
        book(_p_cosDeltaPhi_dy[cat],        5 + 4 * offset, 1, 1);
        book(_p_cosDeltaPhi_pTbar[cat],     6 + 4 * offset, 1, 1);
        book(_p_cosTwoDeltaPhi_dy[cat],    37 + 2 * offset, 1, 1);
        book(_p_cosTwoDeltaPhi_pTbar[cat], 38 + 2 * offset, 1, 1);

        // Gap fraction vs. Q0 and cross-section in dy slices
        _s_gapFrac_Q0.resize(_dy_max);
        for (size_t dyLow = 0; dyLow < _dy_max; ++dyLow ) {
          const string hname("_" + cat + "_dySlice_" + toString(dyLow) + "_" + toString(dyLow+1) + "_Q0");
          { Histo1DPtr tmp; _aux_Q0_dySlices[cat].add(dyLow, dyLow+1, book(tmp, hname, refData(29+dyLow, 1, 1))); }
          { Histo1DPtr tmp; _h_dphi_dySlices[cat].add(dyLow, dyLow+1, book(tmp, 13+(_dy_max*offset)+dyLow, 1, 1)); }
          if (!offset)  book(_s_gapFrac_Q0[dyLow], 29 + dyLow, 1, 1); // only book once
        }

      }

      // Number of jets in rapidity interval
      book(_s_gapFrac_dy,     1, 1, 1);
      book(_s_gapFrac_pTbar,  2, 1, 1);
      book(_p_nGapJets_dy,    3, 1, 1);
      book(_p_nGapJets_pTbar, 4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      for (size_t reg = 0; reg < 2; ++reg) {

        // Retrieve all anti-kt R=0.6 jets 
        const double maxRap = reg? 2.4 : 4.4;
        const Jets& akt6Jets = apply<JetAlg>(event, "AntiKt6JetsWithInvisibles").jetsByPt(Cuts::absrap < maxRap);
        // If there are fewer than 2 jets then bail
        if ( akt6Jets.size() < 2 ) { vetoEvent; }

        // Require jets to be above {60, 50} GeV
        if ( akt6Jets[0].pT() < 60*GeV || akt6Jets[1].pT() < 50*GeV ) { vetoEvent; }

        // Identify gap boundaries
        const double yMin = std::min(akt6Jets[0].rap(), akt6Jets[1].rap());
        const double yMax = std::max(akt6Jets[0].rap(), akt6Jets[1].rap());

        // Determine azimuthal decorrelation quantities
        const double dy = yMax - yMin;
        const double dphi = mapAngle0ToPi(deltaPhi(akt6Jets[0], akt6Jets[1]));
        const double pTbar = 0.5*(akt6Jets[0].pT() + akt6Jets[1].pT())/GeV;

        // Impose minimum dy for the 2011 phase space
        if ( _years[reg] == 2011 && dy < 1.0 ) { vetoEvent; }

        // Determine gap quantities
        size_t nGapJets = 0;
        double maxGapQ0 = 0.;
        const double vetoScale = reg? 30*GeV : 20*GeV;
        for (const Jet& jet : akt6Jets) {
          if (!inRange(jet.rap(), yMin, yMax, OPEN, OPEN))  continue;
          const double pT = jet.pT()/GeV;
          if (pT > vetoScale) { ++nGapJets; }
          if (pT > maxGapQ0) { maxGapQ0 = pT; }
        }

        // Fill histograms
        fillHists(_years[reg], nGapJets, { dy, pTbar, dphi, maxGapQ0 });
      }
      return;
    }

    void fillHists(const size_t region, const size_t nGapJets, const vector<double>& vars) {
      assert(vars.size() == 4);
      const double dy = vars[0];
      const double pTbar = vars[1];
      const double dphi = vars[2];
      const double maxGapQ0 = vars[3];
      // Determine gap category
      vector<string> categories = {"inclusive"};
      if (!nGapJets) { categories += string("gap"); }

      // Fill histograms relevant for comparison with 2010 data
      if (region == _years[0]) {
        // Fill inclusive and gap histograms
        for (const string& cat : categories) {
          _aux_dy[cat]->fill(dy);
          _h_dphi_dySlices[cat].fill(dy, dphi/M_PI);
          _p_cosDeltaPhi_dy[cat]->fill(dy, cos(M_PI - dphi));
          _p_cosTwoDeltaPhi_dy[cat]->fill(dy, cos(2*dphi));
        }
        // Fill profiled nGapJets
        _p_nGapJets_dy->fill(dy, nGapJets);
        // Fill Q0 histograms - can fill multiple points per event
        for (const double Q0 :  _Qnoughts) {
          _aux_Q0_dySlices["inclusive"].fill(dy, Q0);
          if (maxGapQ0 <= Q0) { _aux_Q0_dySlices["gap"].fill(dy, Q0); }
        }

      // Fill histograms relevant for comparison with 2011 data
      } 
      else if (region == _years[1]) {
        // Fill inclusive and gap histograms
        for (const string& cat : categories) {
          _aux_pTbar[cat]->fill(pTbar);
          _p_cosDeltaPhi_pTbar[cat]->fill(pTbar, cos(M_PI - dphi));
          _p_cosTwoDeltaPhi_pTbar[cat]->fill(pTbar, cos(2*dphi));
        }
        // Fill profiled nGapJets
        _p_nGapJets_pTbar->fill(pTbar, nGapJets);
      }
      return;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalise cross-section plots to correct cross-section
      const double ySpan = 1.0; // all dy spans are 1
      const double sf = crossSection() / picobarn / sumOfWeights();
      for (const string& cat : { "inclusive", "gap" }) {
        _h_dphi_dySlices[cat].scale(sf/ySpan/M_PI, this);
        // Create C2/C1 scatter from profiles
        divide(_p_cosTwoDeltaPhi_dy[cat],    _p_cosDeltaPhi_dy[cat],    _h_C2C1_dy[cat]);
        divide(_p_cosTwoDeltaPhi_pTbar[cat], _p_cosDeltaPhi_pTbar[cat], _h_C2C1_pTbar[cat]);
      }

      // Fill simple gap fractions
      efficiency(_aux_dy["gap"], _aux_dy["inclusive"], _s_gapFrac_dy);
      efficiency(_aux_pTbar["gap"], _aux_pTbar["inclusive"], _s_gapFrac_pTbar);

      // Register and fill Q0 gap fractions
      for (size_t dyLow = 0; dyLow < _dy_max; ++dyLow) {
        efficiency(_aux_Q0_dySlices["gap"].histos()[dyLow], _aux_Q0_dySlices["inclusive"].histos()[dyLow], _s_gapFrac_Q0[dyLow]);
      }
    }


  private:

    /// Member variables
    vector<size_t> _years;
    vector<double> _Qnoughts;

    size_t _dy_max;

    /// auxiliary histograms for gap fractions
    map<string, Histo1DPtr> _aux_dy;
    map<string, Histo1DPtr> _aux_pTbar;
    map<string, BinnedHistogram> _aux_Q0_dySlices;

    // Gap fractions
    Scatter2DPtr _s_gapFrac_dy, _s_gapFrac_pTbar;
    vector<Scatter2DPtr> _s_gapFrac_Q0;

    // Number of jets in rapidity interval
    Profile1DPtr _p_nGapJets_dy;
    Profile1DPtr _p_nGapJets_pTbar;

    // Azimuthal moment histograms
    map<string, Profile1DPtr> _p_cosDeltaPhi_dy;
    map<string, Profile1DPtr> _p_cosDeltaPhi_pTbar;
    map<string, Profile1DPtr> _p_cosTwoDeltaPhi_dy;
    map<string, Profile1DPtr> _p_cosTwoDeltaPhi_pTbar;
    map<string, Scatter2DPtr> _h_C2C1_dy;
    map<string, Scatter2DPtr> _h_C2C1_pTbar;

    // Cross-section vs. deltaPhi in deltaY slices
    map<string, BinnedHistogram> _h_dphi_dySlices;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1307243);

}
