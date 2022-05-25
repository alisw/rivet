// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief CMS Z pT and rapidity in Drell-Yan events at 7 TeV
  /// @author Justin Hugon, Luca Perrozzi
  class CMS_2012_I941555 : public Analysis {
  public:

    /// Constructor
    CMS_2012_I941555()
      : Analysis("CMS_2012_I941555")
    {
    }


    /// @name Analysis methods
    //@{

    void init() {

      // Set up projections
      /// @todo Really?: ZFinder zfinder_dressed_mu_pt(-2.1, 2.1, 20, PID::MUON, 60*GeV, 120*GeV, 0.2, false, true);
      FinalState fs;
      Cut cuts = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;
      ZFinder zfinder_dressed_mu_pt(fs, cuts, PID::MUON, 60*GeV, 120*GeV, 0.2);
      declare(zfinder_dressed_mu_pt, "ZFinder_dressed_mu_pt");
      ZFinder zfinder_dressed_el_pt(fs, cuts, PID::ELECTRON, 60*GeV, 120*GeV, 0.1);
      declare(zfinder_dressed_el_pt, "ZFinder_dressed_el_pt");

      ZFinder zfinder_dressed_mu_rap(fs, Cuts::open(), PID::MUON, 60*GeV, 120*GeV, 0.1);
      declare(zfinder_dressed_mu_rap, "ZFinder_dressed_mu_rap");
      ZFinder zfinder_dressed_el_rap(fs, Cuts::open(), PID::ELECTRON, 60*GeV, 120*GeV, 0.1);
      declare(zfinder_dressed_el_rap, "ZFinder_dressed_el_rap");

      // Book histograms
      book(_hist_zrap_mu_dressed      ,1, 1, 1);  // muon "dressed" rapidity
      book(_hist_zrap_el_dressed      ,1, 1, 2);  // electron "dressed" rapidity
      book(_hist_zrap_comb_dressed    ,1, 1, 3);  // electron "dressed" rapidity

      book(_hist_zpt_mu_dressed       ,2, 1, 1);  // muon "dressed" pt
      book(_hist_zpt_el_dressed       ,2, 1, 2);  // electron "dressed" pt
      book(_hist_zpt_comb_dressed     ,2, 1, 3);  // electron "dressed" pt

      book(_hist_zptpeak_mu_dressed   ,3, 1, 1);  // muon "dressed" pt peak
      book(_hist_zptpeak_el_dressed   ,3, 1, 2);  // electron "dressed" pt peak
      book(_hist_zptpeak_comb_dressed ,3, 1, 3);  // electron "dressed" pt peak
    }


    /// Do the analysis
    void analyze(const Event& evt) {
      const double weight = 1.0;

      const ZFinder& zfinder_dressed_mu_rap = apply<ZFinder>(evt, "ZFinder_dressed_mu_rap");
      if (!zfinder_dressed_mu_rap.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_mu_rap.bosons()[0].momentum();
        _hist_zrap_mu_dressed->fill(pZ.rapidity()/GeV, weight);
        _hist_zrap_comb_dressed->fill(pZ.rapidity()/GeV, weight);
      }

      const ZFinder& zfinder_dressed_mu_pt = apply<ZFinder>(evt, "ZFinder_dressed_mu_pt");
      if (!zfinder_dressed_mu_pt.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_mu_pt.bosons()[0].momentum();
        _hist_zpt_mu_dressed->fill(pZ.pT()/GeV, weight);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV, weight);
        if (pZ.pT() < 30*GeV) {
          _hist_zptpeak_mu_dressed->fill(pZ.pT()/GeV, weight);
          _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV, weight);
        }
      }

      const ZFinder& zfinder_dressed_el_rap = apply<ZFinder>(evt, "ZFinder_dressed_el_rap");
      if (!zfinder_dressed_el_rap.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_el_rap.bosons()[0].momentum();
        _hist_zrap_el_dressed->fill(pZ.rapidity()/GeV, weight);
        _hist_zrap_comb_dressed->fill(pZ.rapidity()/GeV, weight);
      }

      const ZFinder& zfinder_dressed_el_pt = apply<ZFinder>(evt, "ZFinder_dressed_el_pt");
      if (!zfinder_dressed_el_pt.bosons().empty()) {
        const FourMomentum pZ = zfinder_dressed_el_pt.bosons()[0].momentum();
        _hist_zpt_el_dressed->fill(pZ.pT()/GeV, weight);
        _hist_zpt_comb_dressed->fill(pZ.pT()/GeV, weight);
        if (pZ.pT() < 30*GeV) {
          _hist_zptpeak_el_dressed->fill(pZ.pT()/GeV, weight);
          _hist_zptpeak_comb_dressed->fill(pZ.pT()/GeV, weight);
        }
      }

    }


    void finalize() {
      normalize(_hist_zrap_mu_dressed);
      normalize(_hist_zpt_mu_dressed);
      normalize(_hist_zptpeak_mu_dressed);

      normalize(_hist_zrap_el_dressed);
      normalize(_hist_zpt_el_dressed);
      normalize(_hist_zptpeak_el_dressed);

      normalize(_hist_zrap_comb_dressed);
      normalize(_hist_zpt_comb_dressed);
      normalize(_hist_zptpeak_comb_dressed);
    }

    //@}


  private:

    Histo1DPtr _hist_zrap_mu_dressed;
    Histo1DPtr _hist_zpt_mu_dressed;
    Histo1DPtr _hist_zptpeak_mu_dressed;

    Histo1DPtr _hist_zrap_el_dressed;
    Histo1DPtr _hist_zpt_el_dressed;
    Histo1DPtr _hist_zptpeak_el_dressed;

    Histo1DPtr _hist_zrap_comb_dressed;
    Histo1DPtr _hist_zpt_comb_dressed;
    Histo1DPtr _hist_zptpeak_comb_dressed;

  };


  RIVET_DECLARE_PLUGIN(CMS_2012_I941555);

}
