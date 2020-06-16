// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief D0 Run II measurement of W charge asymmetry
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2008_S7837160 : public Analysis {
  public:

    /// Default constructor.
    D0_2008_S7837160()
      : Analysis("D0_2008_S7837160")
    {    }


    /// @name Analysis methods
    //@{

    // Book histograms and set up projections
    void init() {
      // Projections
      FinalState fs;
      /// @todo Use separate pT and ETmiss cuts in WFinder
      const WFinder wfe(fs, Cuts::abseta < 5 && Cuts::pT > 25*GeV, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wfe, "WFe");

      // Histograms (temporary +- charge histos and scatters to store the calculated asymmetries)
      for (size_t pmindex = 0; pmindex < 2; ++pmindex) {
        const string suffix = (pmindex == 0) ? "plus" : "minus";
        book(_hs_dsigpm_deta_25_35[pmindex] ,"TMP/dsigpm_deta_25_35_" + suffix, refData(1, 1, 1));
        book(_hs_dsigpm_deta_35[pmindex] ,"TMP/dsigpm_deta_35_" + suffix, refData(1, 1, 2));
        book(_hs_dsigpm_deta_25[pmindex] ,"TMP/dsigpm_deta_25_" + suffix, refData(1, 1, 3));
      }
      book(_h_asym1, 1, 1, 1);
      book(_h_asym2, 1, 1, 2);
      book(_h_asym3, 1, 1, 3);
    }


    /// Do the analysis
    void analyze(const Event & event) {
      const WFinder& wf = apply<WFinder>(event, "WFe");
      if (wf.bosons().size() == 0) {
        MSG_DEBUG("No W candidates found: vetoing");
        vetoEvent;
      }

      // Get the e+- momentum, and an effective charge including the eta sign
      /// @todo Is it correct to multiply the eta sign into the charge to "fold" the plot?
      const FourMomentum p_e = wf.constituentLeptons()[0].momentum();
      const int chg_e = sign(p_e.eta()) * sign(charge(wf.constituentLeptons()[0]));
      assert(chg_e == 1 || chg_e == -1);
      MSG_TRACE("Charged lepton sign = " << chg_e);

      // Fill histos with appropriate +- indexing
      const size_t pmindex = (chg_e > 0) ? 0 : 1;
      if (p_e.Et() < 35*GeV) _hs_dsigpm_deta_25_35[pmindex]->fill(fabs(p_e.eta()));
      else _hs_dsigpm_deta_35[pmindex]->fill(fabs(p_e.eta()));
      _hs_dsigpm_deta_25[pmindex]->fill(fabs(p_e.eta()));
    }


    /// @name Helper functions for constructing asymmetry histograms in finalize()
    //@{
    void calc_asymm(const Histo1DPtr plus, const Histo1DPtr minus, Scatter2DPtr target) {
      divide(*plus - *minus, *plus + *minus, target);
    }
    void calc_asymm(const Histo1DPtr histos[2], Scatter2DPtr target) {
      calc_asymm(histos[0], histos[1], target);
    }
    //@}


    /// @brief Finalize
    ///
    /// Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each ET region
    void finalize() {
      calc_asymm(_hs_dsigpm_deta_25_35, _h_asym1);
      calc_asymm(_hs_dsigpm_deta_35, _h_asym2);
      calc_asymm(_hs_dsigpm_deta_25, _h_asym3);
      _h_asym1->scale(1.,100.);
      _h_asym2->scale(1.,100.);
      _h_asym3->scale(1.,100.);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hs_dsigpm_deta_25_35[2], _hs_dsigpm_deta_35[2], _hs_dsigpm_deta_25[2];
    Scatter2DPtr _h_asym1, _h_asym2, _h_asym3;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2008_S7837160);

}
