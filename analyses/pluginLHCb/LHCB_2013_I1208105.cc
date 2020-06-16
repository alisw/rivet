// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class LHCB_2013_I1208105 : public Analysis {
  public:

    LHCB_2013_I1208105()
      : Analysis("LHCB_2013_I1208105")
    {   }


    void init() {
      // Projections
      declare(FinalState((Cuts::etaIn(1.9, 4.9))), "forwardFS");
      declare(FinalState((Cuts::etaIn(-3.5,-1.5))), "backwardFS");
      declare(ChargedFinalState((Cuts::etaIn(1.9, 4.9))), "forwardCFS");
      declare(ChargedFinalState((Cuts::etaIn(-3.5,-1.5))), "backwardCFS");

      // Histos
      book(_s_chEF_minbias, 1, 1, 1, true);
      book(_s_chEF_hard, 2, 1, 1, true);
      book(_s_chEF_diff, 3, 1, 1, true);
      book(_s_chEF_nondiff, 4, 1, 1, true);
      book(_s_totEF_minbias, 5, 1, 1, true);
      book(_s_totEF_hard, 6, 1, 1, true);
      book(_s_totEF_diff, 7, 1, 1, true);
      book(_s_totEF_nondiff, 8, 1, 1, true);

      // Temporary profiles and histos
      /// @todo Convert to declared/registered temp histos
      book(_tp_chEF_minbias, "TMP/chEF_minbias", refData(1,1,1));
      book(_tp_chEF_hard, "TMP/chEF_hard", refData(2,1,1));
      book(_tp_chEF_diff, "TMP/chEF_diff", refData(3,1,1));
      book(_tp_chEF_nondiff, "TMP/chEF_nondiff", refData(4,1,1));
      book(_tp_totEF_minbias, "TMP/totEF_minbias", refData(5,1,1));
      book(_tp_totEF_hard, "TMP/totEF_hard", refData(6,1,1));
      book(_tp_totEF_diff, "TMP/totEF_diff", refData(7,1,1));
      book(_tp_totEF_nondiff, "TMP/totEF_nondiff", refData(8,1,1));

      book(_th_chN_minbias, "TMP/chN_minbias", refData(1,1,1));
      book(_th_chN_hard, "TMP/chN_hard", refData(2,1,1));
      book(_th_chN_diff, "TMP/chN_diff", refData(3,1,1));
      book(_th_chN_nondiff, "TMP/chN_nondiff", refData(4,1,1));
      book(_th_totN_minbias, "TMP/totN_minbias", refData(5,1,1));
      book(_th_totN_hard, "TMP/totN_hard", refData(6,1,1));
      book(_th_totN_diff, "TMP/totN_diff", refData(7,1,1));
      book(_th_totN_nondiff, "TMP/totN_nondiff", refData(8,1,1));

      // Counters
      book(_mbSumW, "TMP/mbSumW");
      book(_hdSumW, "TMP/hdSumW");
      book(_dfSumW, "TMP/dfSumW");
      book(_ndSumW, "TMP/ndSumW");
      book(_mbchSumW, "TMP/mbchSumW");
      book(_hdchSumW, "TMP/hdchSumW");
      book(_dfchSumW, "TMP/dfchSumW");
      book(_ndchSumW, "TMP/ndchSumW");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ffs = apply<FinalState>(event, "forwardFS");
      const FinalState& bfs = apply<FinalState>(event, "backwardFS");
      const ChargedFinalState& fcfs = apply<ChargedFinalState>(event, "forwardCFS");
      const ChargedFinalState& bcfs = apply<ChargedFinalState>(event, "backwardCFS");

      // Veto this event completely if there are no forward *charged* particles
      if (fcfs.empty()) vetoEvent;

      // Charged and neutral version
      {
        // Decide empirically if this is a "hard" or "diffractive" event
        bool ishardEvt = false;
        for (const Particle& p : ffs.particles()) {
          if (p.pT() > 3.0*GeV) { ishardEvt = true; break; }
        }
        // Decide empirically if this is a "diffractive" event
        /// @todo Can be "diffractive" *and* "hard"?
        bool isdiffEvt = (bfs.size() == 0);

        // Update event-type weight counters
        _mbSumW->fill();
        (isdiffEvt ? _dfSumW : _ndSumW)->fill();
        if (ishardEvt) _hdSumW->fill();

        // Plot energy flow
        for (const Particle& p : ffs.particles()) {
          const double eta = p.eta();
          const double energy = p.E();
          _tp_totEF_minbias->fill(eta, energy);
          _th_totN_minbias->fill(eta);
          if (ishardEvt) {
            _tp_totEF_hard->fill(eta, energy);
            _th_totN_hard->fill(eta);
          }
          if (isdiffEvt) {
            _tp_totEF_diff->fill(eta, energy);
            _th_totN_diff->fill(eta);
          } else {
            _tp_totEF_nondiff->fill(eta, energy);
            _th_totN_nondiff->fill(eta);
          }
        }
      }


      // Charged-only version
      {
        bool ishardEvt = false;
        for (const Particle& p : fcfs.particles()) {
          if (p.pT() > 3.0*GeV) { ishardEvt = true; break; }
        }
        // Decide empirically if this is a "diffractive" event
        /// @todo Can be "diffractive" *and* "hard"?
        bool isdiffEvt = (bcfs.size() == 0);

        // Update event-type weight counters
        _mbchSumW->fill();
        (isdiffEvt ? _dfchSumW : _ndchSumW)->fill();
        if (ishardEvt) _hdchSumW->fill();

        // Plot energy flow
        for (const Particle& p : fcfs.particles()) {
          const double eta = p.eta();
          const double energy = p.E();
          _tp_chEF_minbias->fill(eta, energy);
          _th_chN_minbias->fill(eta);
          if (ishardEvt) {
            _tp_chEF_hard->fill(eta, energy);
            _th_chN_hard->fill(eta);
          }
          if (isdiffEvt) {
            _tp_chEF_diff->fill(eta, energy);
            _th_chN_diff->fill(eta);
          } else {
            _tp_chEF_nondiff->fill(eta, energy);
            _th_chN_nondiff->fill(eta);
          }
        }
      }

    }


    void finalize() {
      for (size_t i = 0; i < _s_totEF_minbias->numPoints(); ++i) {
        const double val = _tp_totEF_minbias->bin(i).mean() * _th_totN_minbias->bin(i).height();
        const double err = (_tp_totEF_minbias->bin(i).mean() * _th_totN_minbias->bin(i).heightErr() +
                            _tp_totEF_minbias->bin(i).stdErr() * _th_totN_minbias->bin(i).height());
        _s_totEF_minbias->point(i).setY(val/_mbSumW->val(), err/_mbSumW->val());
      }
      for (size_t i = 0; i < _s_totEF_hard->numPoints(); ++i) {
        const double val = _tp_totEF_hard->bin(i).mean() * _th_totN_hard->bin(i).height();
        const double err = (_tp_totEF_hard->bin(i).mean() * _th_totN_hard->bin(i).heightErr() +
                            _tp_totEF_hard->bin(i).stdErr() * _th_totN_hard->bin(i).height());
        _s_totEF_hard->point(i).setY(val/_hdSumW->val(), err/_hdSumW->val());
      }
      for (size_t i = 0; i < _s_totEF_diff->numPoints(); ++i) {
        const double val = _tp_totEF_diff->bin(i).mean() * _th_totN_diff->bin(i).height();
        const double err = (_tp_totEF_diff->bin(i).mean() * _th_totN_diff->bin(i).heightErr() +
                                   _tp_totEF_diff->bin(i).stdErr() * _th_totN_diff->bin(i).height());
        _s_totEF_diff->point(i).setY(val/_dfSumW->val(), err/_dfSumW->val());
      }
      for (size_t i = 0; i < _s_totEF_nondiff->numPoints(); ++i) {
        const double val = _tp_totEF_nondiff->bin(i).mean() * _th_totN_nondiff->bin(i).height();
        const double err = (_tp_totEF_nondiff->bin(i).mean() * _th_totN_nondiff->bin(i).heightErr() +
                            _tp_totEF_nondiff->bin(i).stdErr() * _th_totN_nondiff->bin(i).height());
        _s_totEF_nondiff->point(i).setY(val/_ndSumW->val(), err/_ndSumW->val());
      }
      for (size_t i = 0; i < _s_chEF_minbias->numPoints(); ++i) {
        const double val = _tp_chEF_minbias->bin(i).mean() * _th_chN_minbias->bin(i).height();
        const double err = (_tp_chEF_minbias->bin(i).mean() * _th_chN_minbias->bin(i).heightErr() +
                            _tp_chEF_minbias->bin(i).stdErr() * _th_chN_minbias->bin(i).height());
        _s_chEF_minbias->point(i).setY(val/_mbchSumW->val(), err/_mbchSumW->val());
      }
      for (size_t i = 0; i < _s_chEF_hard->numPoints(); ++i) {
        const double val = _tp_chEF_hard->bin(i).mean() * _th_chN_hard->bin(i).height();
        const double err = (_tp_chEF_hard->bin(i).mean() * _th_chN_hard->bin(i).heightErr() +
                            _tp_chEF_hard->bin(i).stdErr() * _th_chN_hard->bin(i).height());
        _s_chEF_hard->point(i).setY(val/_hdchSumW->val(), err/_hdchSumW->val());
      }
      for (size_t i = 0; i < _s_chEF_diff->numPoints(); ++i) {
        const double val = _tp_chEF_diff->bin(i).mean() * _th_chN_diff->bin(i).height();
        const double err = (_tp_chEF_diff->bin(i).mean() * _th_chN_diff->bin(i).heightErr() +
                            _tp_chEF_diff->bin(i).stdErr() * _th_chN_diff->bin(i).height());
        _s_chEF_diff->point(i).setY(val/_dfchSumW->val(), err/_dfchSumW->val());
      }
      for (size_t i = 0; i < _s_chEF_nondiff->numPoints(); ++i) {
        const double val = _tp_chEF_nondiff->bin(i).mean() * _th_chN_nondiff->bin(i).height();
        const double err = (_tp_chEF_nondiff->bin(i).mean() * _th_chN_nondiff->bin(i).heightErr() +
                            _tp_chEF_nondiff->bin(i).stdErr() * _th_chN_nondiff->bin(i).height());
        _s_chEF_nondiff->point(i).setY(val/_ndchSumW->val(), err/_ndchSumW->val());
      }
    }


  private:

    /// @name Histograms and counters
    ///
    /// @note Histograms correspond to charged and total EF for each class of events:
    ///  minimum bias, hard scattering, diffractive enriched and non-diffractive enriched.
    //@{

    // Scatters to be filled in finalize with 1/d_eta <N(eta)><E(eta)>
    Scatter2DPtr _s_totEF_minbias, _s_totEF_hard, _s_totEF_diff, _s_totEF_nondiff;
    Scatter2DPtr _s_chEF_minbias, _s_chEF_hard, _s_chEF_diff, _s_chEF_nondiff;

    // Temp profiles containing <E(eta)>
    Profile1DPtr _tp_totEF_minbias, _tp_totEF_hard, _tp_totEF_diff, _tp_totEF_nondiff;
    Profile1DPtr _tp_chEF_minbias, _tp_chEF_hard, _tp_chEF_diff, _tp_chEF_nondiff;

    // Temp profiles containing <N(eta)>
    Histo1DPtr _th_totN_minbias, _th_totN_hard, _th_totN_diff, _th_totN_nondiff;
    Histo1DPtr _th_chN_minbias, _th_chN_hard, _th_chN_diff, _th_chN_nondiff;

    // Sums of weights (~ #events) in each event class
    CounterPtr _mbSumW, _hdSumW, _dfSumW, _ndSumW;
    CounterPtr _mbchSumW, _hdchSumW, _dfchSumW, _ndchSumW;

    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2013_I1208105);

}
