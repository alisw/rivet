// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2017_I1598460 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1598460);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs, "FinalState");
      FastJets fj07(fs, FastJets::ANTIKT, 0.7);
      declare(fj07, "Jets");
      /// @todo Book histograms here, e.g.:
      for (int i = 0; i < 6; i++) {
        Histo1DPtr tmp; _h_ybys.push_back(book(tmp, 2*i+1,1,1));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FastJets& fj = apply<FastJets>(event,"Jets");
      const Jets& jets = fj.jetsByPt(Cuts::pt > 50*GeV && Cuts::absrap < 5);

      // Require two jets
      if (jets.size() < 2) vetoEvent;

      // Veto events if one of two leading jets |y|>3.0, otherwise
      // the subleading jets can become the leading jets through jet selection
      if (jets[0].absrap() > 3 || jets[1].absrap() > 3) vetoEvent;

      double ystar = 0.5 * std::abs(jets[0].rap() - jets[1].rap());
      double yboost = 0.5 * std::abs(jets[0].rap() + jets[1].rap());
      double ptavg = 0.5 * (jets[0].pT() + jets[1].pT());

      // Compute index of histogram to be filled: yb0ys0 --> 0 yb0ys1 --> 1 ...
      size_t i = (size_t)yboost;
      size_t j = (size_t)ystar;
      size_t idx = j + 3*i - i*(i-1)/2;

      _h_ybys[idx]->fill(ptavg/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_ybys, crossSection()/sumOfWeights());
    }

    //@}


    vector<Histo1DPtr> _h_ybys;

  };


  DECLARE_RIVET_PLUGIN(CMS_2017_I1598460);

}
