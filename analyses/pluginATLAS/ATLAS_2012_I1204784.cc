// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  


  /// ATLAS Z phi* measurement
  class ATLAS_2012_I1204784 : public Analysis {
    public:

      /// Constructor
      ATLAS_2012_I1204784()
        : Analysis("ATLAS_2012_I1204784")
      {      }


    public:

      /// Book histograms and initialise projections before the run
      void init() {
        FinalState fs;
        Cut cuts = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
        ZFinder zfinder_dressed_el(fs, cuts, PID::ELECTRON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
        declare(zfinder_dressed_el, "ZFinder_dressed_el");
        ZFinder zfinder_bare_el(fs, cuts, PID::ELECTRON, 66*GeV, 116*GeV, 0.0, ZFinder::ClusterPhotons::NONE);
        declare(zfinder_bare_el, "ZFinder_bare_el");
        ZFinder zfinder_dressed_mu(fs, cuts, PID::MUON, 66*GeV, 116*GeV, 0.1, ZFinder::ClusterPhotons::NODECAY);
        declare(zfinder_dressed_mu, "ZFinder_dressed_mu");
        ZFinder zfinder_bare_mu(fs, cuts, PID::MUON, 66*GeV, 116*GeV, 0.0, ZFinder::ClusterPhotons::NONE);
        declare(zfinder_bare_mu, "ZFinder_bare_mu");

        // Book histograms
        // Single-differential plots
        book(_hist_zphistar_el_bare ,1, 1, 1);
        book(_hist_zphistar_mu_bare ,1, 1, 2);
        book(_hist_zphistar_el_dressed ,2, 1, 1);
        book(_hist_zphistar_mu_dressed ,2, 1, 2);

        // Double-differential plots
        {Histo1DPtr tmp; _h_phistar_el_bare.add(0.0, 0.8, book(tmp, 3, 1, 1));}
        {Histo1DPtr tmp; _h_phistar_el_bare.add(0.8, 1.6, book(tmp, 3, 1, 2));}
        {Histo1DPtr tmp; _h_phistar_el_bare.add(1.6, 10.0, book(tmp, 3, 1, 3));}

        {Histo1DPtr tmp; _h_phistar_el_dressed.add(0.0, 0.8, book(tmp, 3, 2, 1));}
        {Histo1DPtr tmp; _h_phistar_el_dressed.add(0.8, 1.6, book(tmp, 3, 2, 2));}
        {Histo1DPtr tmp; _h_phistar_el_dressed.add(1.6, 10.0, book(tmp, 3, 2, 3));}

        {Histo1DPtr tmp; _h_phistar_mu_bare.add(0.0, 0.8, book(tmp, 4, 1, 1));}
        {Histo1DPtr tmp; _h_phistar_mu_bare.add(0.8, 1.6, book(tmp, 4, 1, 2));}
        {Histo1DPtr tmp; _h_phistar_mu_bare.add(1.6, 10.0, book(tmp, 4, 1, 3));}

        {Histo1DPtr tmp; _h_phistar_mu_dressed.add(0.0, 0.8, book(tmp, 4, 2, 1));}
        {Histo1DPtr tmp; _h_phistar_mu_dressed.add(0.8, 1.6, book(tmp, 4, 2, 2));}
        {Histo1DPtr tmp; _h_phistar_mu_dressed.add(1.6, 10.0, book(tmp, 4, 2, 3));}
      }


      /// Perform the per-event analysis
      void analyze(const Event& event) {
        const double weight = 1.0;

        const ZFinder& zfinder_dressed_el = apply<ZFinder>(event, "ZFinder_dressed_el");
        const ZFinder& zfinder_bare_el = apply<ZFinder>(event, "ZFinder_bare_el");
        const ZFinder& zfinder_dressed_mu = apply<ZFinder>(event, "ZFinder_dressed_mu");
        const ZFinder& zfinder_bare_mu = apply<ZFinder>(event, "ZFinder_bare_mu");

        fillPlots(zfinder_dressed_el, _hist_zphistar_el_dressed, _h_phistar_el_dressed, weight);
        fillPlots(zfinder_bare_el, _hist_zphistar_el_bare, _h_phistar_el_bare, weight);
        fillPlots(zfinder_dressed_mu, _hist_zphistar_mu_dressed, _h_phistar_mu_dressed, weight);
        fillPlots(zfinder_bare_mu, _hist_zphistar_mu_bare, _h_phistar_mu_bare, weight);
      }


      void fillPlots(const ZFinder& zfind, Histo1DPtr hist, BinnedHistogram& binnedHist, double weight) {
        if (zfind.bosons().size() != 1) return;
        Particles leptons = sortBy(zfind.constituents(), cmpMomByPt);

        const FourMomentum lminus = leptons[0].charge() < 0 ? leptons[0].momentum() : leptons[1].momentum();
        const FourMomentum lplus = leptons[0].charge() < 0 ? leptons[1].momentum() : leptons[0].momentum();

        const double phi_acop = M_PI - deltaPhi(lminus, lplus);
        const double costhetastar = tanh((lminus.eta()-lplus.eta())/2.0);
        const double sin2thetastar = (costhetastar <= 1) ? 1.0 - sqr(costhetastar) : 0;
        const double phistar = tan(phi_acop/2.0) * sqrt(sin2thetastar);
        hist->fill(phistar, weight);

        binnedHist.fill(zfind.bosons()[0].absrap(), phistar, weight);
      }


      /// Normalise histograms etc., after the run
      void finalize() {
        normalize(_hist_zphistar_el_dressed);
        normalize(_hist_zphistar_el_bare);
        normalize(_hist_zphistar_mu_dressed);
        normalize(_hist_zphistar_mu_bare);

        for (Histo1DPtr hist : _h_phistar_mu_dressed.histos()) { normalize(hist); }
        for (Histo1DPtr hist : _h_phistar_mu_bare.histos()) { normalize(hist); }
        for (Histo1DPtr hist : _h_phistar_el_bare.histos()) { normalize(hist); }
        for (Histo1DPtr hist : _h_phistar_el_dressed.histos()) { normalize(hist); }
      }

      //@}


    private:

      BinnedHistogram _h_phistar_mu_dressed;
      BinnedHistogram _h_phistar_mu_bare;
      BinnedHistogram _h_phistar_el_dressed;
      BinnedHistogram _h_phistar_el_bare;

      Histo1DPtr _hist_zphistar_el_dressed;
      Histo1DPtr _hist_zphistar_el_bare;

      Histo1DPtr _hist_zphistar_mu_dressed;
      Histo1DPtr _hist_zphistar_mu_bare;

      Histo1DPtr _hist_zphistar_el_bare_1;
      Histo1DPtr _hist_zphistar_el_bare_2;
      Histo1DPtr _hist_zphistar_el_bare_3;

      Histo1DPtr _hist_zphistar_el_dressed_1;
      Histo1DPtr _hist_zphistar_el_dressed_2;
      Histo1DPtr _hist_zphistar_el_dressed_3;

      Histo1DPtr _hist_zphistar_mu_bare_1;
      Histo1DPtr _hist_zphistar_mu_bare_2;
      Histo1DPtr _hist_zphistar_mu_bare_3;

      Histo1DPtr _hist_zphistar_mu_dressed_1;
      Histo1DPtr _hist_zphistar_mu_dressed_2;
      Histo1DPtr _hist_zphistar_mu_dressed_3;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2012_I1204784);

}
