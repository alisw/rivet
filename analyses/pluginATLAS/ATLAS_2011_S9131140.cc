// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief ATLAS Z pT in Drell-Yan events at 7 TeV
  ///
  /// @author Elena Yatsenko, Judith Katzy
  class ATLAS_2011_S9131140 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_S9131140);


    /// @name Analysis methods
    //@{

    void init() {

      // Set up projections
      FinalState fs;
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;

      ZFinder zfinder_dressed_el(fs, cut, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_dressed_el, "ZFinder_dressed_el");
      ZFinder zfinder_bare_el(fs, cut, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.0, ZFinder::ClusterPhotons::NONE);
      declare(zfinder_bare_el, "ZFinder_bare_el");
      ZFinder zfinder_dressed_mu(fs, cut, PID::MUON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
      declare(zfinder_dressed_mu, "ZFinder_dressed_mu");
      ZFinder zfinder_bare_mu(fs, cut, PID::MUON, 66.0*GeV, 116.0*GeV, 0.0, ZFinder::ClusterPhotons::NONE);
      declare(zfinder_bare_mu, "ZFinder_bare_mu");

      // Book histograms
      book(_hist_zpt_el_dressed     ,1, 1, 2);  // electron "dressed"
      book(_hist_zpt_el_bare        ,1, 1, 3);  // electron "bare"
      book(_hist_zpt_mu_dressed     ,2, 1, 2);  // muon "dressed"
      book(_hist_zpt_mu_bare        ,2, 1, 3);  // muon "bare"

      book(_sumw_el_bare, "_sumw_el_bare");
      book(_sumw_el_dressed, "_sumw_el_dressed");
      book(_sumw_mu_bare, "_sumw_mu_bare");
      book(_sumw_mu_dressed, "_sumw_mu_dressed");
    }


    /// Do the analysis
    void analyze(const Event& evt) {
      const ZFinder& zfinder_dressed_el = apply<ZFinder>(evt, "ZFinder_dressed_el");
      if (!zfinder_dressed_el.bosons().empty()) {
        _sumw_el_dressed->fill();
        const FourMomentum pZ = zfinder_dressed_el.bosons()[0].momentum();
        _hist_zpt_el_dressed->fill(pZ.pT()/GeV);
      }

      const ZFinder& zfinder_bare_el = apply<ZFinder>(evt, "ZFinder_bare_el");
      if (!zfinder_bare_el.bosons().empty()) {
        _sumw_el_bare->fill();
	    const FourMomentum pZ = zfinder_bare_el.bosons()[0].momentum();
        _hist_zpt_el_bare->fill(pZ.pT()/GeV);
      }

      const ZFinder& zfinder_dressed_mu = apply<ZFinder>(evt, "ZFinder_dressed_mu");
      if (!zfinder_dressed_mu.bosons().empty()) {
        _sumw_mu_dressed->fill();
        const FourMomentum pZ = zfinder_dressed_mu.bosons()[0].momentum();
        _hist_zpt_mu_dressed->fill(pZ.pT()/GeV);
      }

      const ZFinder& zfinder_bare_mu = apply<ZFinder>(evt, "ZFinder_bare_mu");
      if (!zfinder_bare_mu.bosons().empty()) {
        _sumw_mu_bare->fill();
        const FourMomentum pZ = zfinder_bare_mu.bosons()[0].momentum();
        _hist_zpt_mu_bare->fill(pZ.pT()/GeV);
      }

    }


    void finalize() {
      if (_sumw_el_dressed->val() != 0) scale(_hist_zpt_el_dressed, 1/ *_sumw_el_dressed);
      if (_sumw_el_bare->val()    != 0) scale(_hist_zpt_el_bare,    1/ *_sumw_el_bare);
      if (_sumw_mu_dressed->val() != 0) scale(_hist_zpt_mu_dressed, 1/ *_sumw_mu_dressed);
      if (_sumw_mu_bare->val()    != 0) scale(_hist_zpt_mu_bare,    1/ *_sumw_mu_bare);
    }

    //@}


    private:

	CounterPtr _sumw_el_bare, _sumw_el_dressed;
	CounterPtr _sumw_mu_bare, _sumw_mu_dressed;

	Histo1DPtr _hist_zpt_el_dressed;
	Histo1DPtr _hist_zpt_el_bare;
	Histo1DPtr _hist_zpt_mu_dressed;
	Histo1DPtr _hist_zpt_mu_bare;
  };



  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_S9131140, ATLAS_2011_I917931);

}
