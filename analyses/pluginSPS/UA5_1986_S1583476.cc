// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/TriggerUA5.hh"

namespace Rivet {


  /// @brief UA5 \f$ \eta \f$ distributions at 200 and 900 GeV
  class UA5_1986_S1583476 : public Analysis {
  public:

    /// Constructor
    UA5_1986_S1583476() : Analysis("UA5_1986_S1583476") {
    }


    /// @name Analysis methods
    //@{

    /// Set up projections and histograms
    void init() {
      declare(TriggerUA5(), "Trigger");
      declare(Beam(), "Beams");
      declare(ChargedFinalState((Cuts::etaIn(-5.0, 5.0))), "CFS50");

      // Histograms
      if (fuzzyEquals(sqrtS()/GeV, 200.0, 1E-4)) {
        book(_hist_eta_nsd       ,1,1,1);
        book(_hist_eta_inelastic ,1,1,2);
        _hists_eta_nsd.resize(6);
        for (int i = 1; i <= 6; ++i) {
          _sumWn.push_back({});
          book(_sumWn.back(), "TMP/sumWn"+to_str(i));
          book(_hists_eta_nsd[i-1],2,1,i);
        }
      } else if (fuzzyEquals(sqrtS()/GeV, 900.0, 1E-4)) {
        book(_hist_eta_nsd       ,1,1,3);
        book(_hist_eta_inelastic ,1,1,4);
        _hists_eta_nsd.resize(9);
        for (int i = 1; i <= 9; ++i) {
          _sumWn.push_back({});
          book(_sumWn.back(), "TMP/sumWn"+to_str(i));
          book(_hists_eta_nsd[i-1],3,1,i);
        }
      }
      book(_sumWTrig, "sumWtrig");
      book(_sumWTrigNSD, "sumWtrigNSD");

    }


    /// Fill eta histograms (in Nch bins)
    void analyze(const Event& event) {
      // Trigger
      const TriggerUA5& trigger = apply<TriggerUA5>(event, "Trigger");
      if (!trigger.sdDecision()) vetoEvent;
      const bool isNSD = trigger.nsdDecision();

      // Get the index corresponding to the max Nch range histo/sum(w) vector index
      const ChargedFinalState& cfs50 = apply<ChargedFinalState>(event, "CFS50");
      const int numP = cfs50.size();
      const int ni = (int)floor(static_cast<float>(numP-2)/10.0);
      const int num_idx = min(ni, (int)_sumWn.size()-1);
      MSG_TRACE("Multiplicity index: " << numP << " charged particles -> #" << num_idx);

      // Update weights
      _sumWTrig->fill();
      if (isNSD) {
        _sumWTrigNSD->fill();
        if (num_idx >= 0) _sumWn[num_idx]->fill();
      }

      // Fill histos
      for (const Particle& p : cfs50.particles()) {
        const double eta = p.abseta();
        _hist_eta_inelastic->fill(eta);
        if (isNSD) {
          _hist_eta_nsd->fill(eta);
          if (num_idx >= 0) _hists_eta_nsd[num_idx]->fill(eta);
        }
      }
    }


    /// Scale histos
    void finalize() {
      MSG_DEBUG("sumW_NSD,inel = " << _sumWTrigNSD->val() << ", " << _sumWTrig->val());
      scale(_hist_eta_nsd, 0.5 / *_sumWTrigNSD);
      scale(_hist_eta_inelastic, 0.5 / *_sumWTrig);
      //
      for (size_t i = 0; i < _hists_eta_nsd.size(); ++i) {
        MSG_DEBUG("sumW[n] = " << _sumWn[i]->val());
        scale(_hists_eta_nsd[i], 0.5 / *_sumWn[i]);
      }
    }


  private:

    /// @name Weight counters
    //@{
    CounterPtr _sumWTrig;
    CounterPtr _sumWTrigNSD;
    vector<CounterPtr> _sumWn;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _hist_eta_nsd;
    Histo1DPtr _hist_eta_inelastic;
    vector<Histo1DPtr> _hists_eta_nsd;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(UA5_1986_S1583476);

}
