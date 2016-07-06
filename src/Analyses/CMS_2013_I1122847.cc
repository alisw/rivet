// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  using namespace Cuts;

  class CMS_2013_I1122847 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1122847()
      : Analysis("CMS_2013_I1122847")  {}

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;

      Cut cuts_mu = etaIn(-2.4, 2.4) & (pT >= 20.0*GeV);
      ZFinder zfinder_mu(fs, cuts_mu, PID::MUON, 40.0*GeV, MAXDOUBLE,
                         0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      addProjection(zfinder_mu, "zfinder_mu");

      Cut cuts_el = (pT >= 20.0*GeV) & ((abseta < 1.447) | ((abseta > 1.57) & (abseta < 2.4)));
      ZFinder zfinder_el(fs, cuts_el, PID::ELECTRON, 40.0*GeV, MAXDOUBLE,
                         0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      addProjection(zfinder_el, "zfinder_el");


      /// histograms
      // dimuon
      _hist_mm_100_num = Histo1D(refData(1, 1, 1));
      _hist_mm_125_num = Histo1D(refData(1, 1, 2));
      _hist_mm_150_num = Histo1D(refData(1, 1, 3));
      _hist_mm_240_num = Histo1D(refData(1, 1, 4));

      _hist_mm_100_den = Histo1D(refData(1, 1, 1));
      _hist_mm_125_den = Histo1D(refData(1, 1, 2));
      _hist_mm_150_den = Histo1D(refData(1, 1, 3));
      _hist_mm_240_den = Histo1D(refData(1, 1, 4));

      // dielectron
      _hist_ee_100_num = Histo1D(refData(2, 1, 1));
      _hist_ee_125_num = Histo1D(refData(2, 1, 2));
      _hist_ee_150_num = Histo1D(refData(2, 1, 3));
      _hist_ee_240_num = Histo1D(refData(2, 1, 4));

      _hist_ee_100_den = Histo1D(refData(2, 1, 1));
      _hist_ee_125_den = Histo1D(refData(2, 1, 2));
      _hist_ee_150_den = Histo1D(refData(2, 1, 3));
      _hist_ee_240_den = Histo1D(refData(2, 1, 4));

      // dilepton
      _hist_ll_100_num = Histo1D(refData(3, 1, 1));
      _hist_ll_125_num = Histo1D(refData(3, 1, 2));
      _hist_ll_150_num = Histo1D(refData(3, 1, 3));
      _hist_ll_240_num = Histo1D(refData(3, 1, 4));

      _hist_ll_100_den = Histo1D(refData(3, 1, 1));
      _hist_ll_125_den = Histo1D(refData(3, 1, 2));
      _hist_ll_150_den = Histo1D(refData(3, 1, 3));
      _hist_ll_240_den = Histo1D(refData(3, 1, 4));
    }

    double cosThetaCS(const Particle& l1, const Particle&  l2) {

      FourMomentum mom1 = l1.momentum();
      FourMomentum mom2 = l2.momentum();

      double Q = FourMomentum(mom1 + mom2).mass();
      double QT = sqrt(pow((mom1.px() + mom2.px()), 2) + pow((mom1.py() + mom2.py()), 2));

      double P1p = 0.0;
      double P1m = 0.0;
      double P2p = 0.0;
      double P2m = 0.0;


      if (l1.pid() > 0) {
        P1p = (mom1.E() + mom1.pz()) / sqrt(2.0);
        P1m = (mom1.E() - mom1.pz()) / sqrt(2.0);
        P2p = (mom2.E() + mom2.pz()) / sqrt(2.0);
        P2m = (mom2.E() - mom2.pz()) / sqrt(2.0);
      } else if (l1.pid() < 0) {
        P1p = (mom2.E() + mom2.pz()) / sqrt(2.0);
        P1m = (mom2.E() - mom2.pz()) / sqrt(2.0);
        P2p = (mom1.E() + mom1.pz()) / sqrt(2.0);
        P2m = (mom1.E() - mom1.pz()) / sqrt(2.0);
      }

      double QZ = mom1.pz() + mom2.pz();
      double cosThetaCS = (2 / (Q * sqrt(Q * Q + QT * QT))) * (P1p * P2m - P1m * P2p);
      if (QZ < 0.0)
        cosThetaCS = -cosThetaCS;

      return cosThetaCS;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const ZFinder& zfinder_el = applyProjection<ZFinder>(event, "zfinder_el");
      if (zfinder_el.bosons().size() > 0) {
        const Particle& z  = zfinder_el.bosons()[0];
        const Particle& l1 = zfinder_el.constituents()[0];
        const Particle& l2 = zfinder_el.constituents()[1];

        // Prepare variables for filling
        double rap = z.momentum().absrap();
        double costhetacs = cosThetaCS(l1, l2);

        double sign = 1.0;
        if (costhetacs < 0.) sign = -1.0;
        // Fill the histograms
        if (rap < 1.0) {
          _hist_ee_100_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_100_num.fill(z.momentum().mass(), weight * sign);
          _hist_ee_100_den.fill(z.momentum().mass(), weight);
          _hist_ll_100_den.fill(z.momentum().mass(), weight);
        } else if (rap < 1.25) {
          _hist_ee_125_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_125_num.fill(z.momentum().mass(), weight * sign);
          _hist_ee_125_den.fill(z.momentum().mass(), weight);
          _hist_ll_125_den.fill(z.momentum().mass(), weight);
        } else if (rap < 1.50) {
          _hist_ee_150_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_150_num.fill(z.momentum().mass(), weight * sign);
          _hist_ee_150_den.fill(z.momentum().mass(), weight);
          _hist_ll_150_den.fill(z.momentum().mass(), weight);
        } else if (rap < 2.40) {
          _hist_ee_240_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_240_num.fill(z.momentum().mass(), weight * sign);
          _hist_ee_240_den.fill(z.momentum().mass(), weight);
          _hist_ll_240_den.fill(z.momentum().mass(), weight);
        }
      }

      const ZFinder& zfinder_mu = applyProjection<ZFinder>(event, "zfinder_mu");
      if (zfinder_mu.bosons().size() > 0) {
        const Particle& z  = zfinder_mu.bosons()[0];
        const Particle& l1 = zfinder_mu.constituents()[0];
        const Particle& l2 = zfinder_mu.constituents()[1];

        // Prepare variables for filling
        double rap = z.momentum().absrap();
        double costhetacs = cosThetaCS(l1, l2);

        double sign = 1.0;
        if (costhetacs < 0.) sign = -1.0;

        // Fill the histograms
        if (rap < 1.0) {
          _hist_mm_100_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_100_num.fill(z.momentum().mass(), weight * sign);
          _hist_mm_100_den.fill(z.momentum().mass(), weight);
          _hist_ll_100_den.fill(z.momentum().mass(), weight);
        } else if (rap < 1.25) {
          _hist_mm_125_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_125_num.fill(z.momentum().mass(), weight * sign);
          _hist_mm_125_den.fill(z.momentum().mass(), weight);
          _hist_ll_125_den.fill(z.momentum().mass(), weight);
        } else if (rap < 1.50) {
          _hist_mm_150_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_150_num.fill(z.momentum().mass(), weight * sign);
          _hist_mm_150_den.fill(z.momentum().mass(), weight);
          _hist_ll_150_den.fill(z.momentum().mass(), weight);
        } else if (rap < 2.40) {
          _hist_mm_240_num.fill(z.momentum().mass(), weight * sign);
          _hist_ll_240_num.fill(z.momentum().mass(), weight * sign);
          _hist_mm_240_den.fill(z.momentum().mass(), weight);
          _hist_ll_240_den.fill(z.momentum().mass(), weight);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      divide(_hist_mm_100_num, _hist_mm_100_den, bookScatter2D(1, 1, 1));
      divide(_hist_mm_125_num, _hist_mm_125_den, bookScatter2D(1, 1, 2));
      divide(_hist_mm_150_num, _hist_mm_150_den, bookScatter2D(1, 1, 3));
      divide(_hist_mm_240_num, _hist_mm_240_den, bookScatter2D(1, 1, 4));

      divide(_hist_ee_100_num, _hist_ee_100_den, bookScatter2D(2, 1, 1));
      divide(_hist_ee_125_num, _hist_ee_125_den, bookScatter2D(2, 1, 2));
      divide(_hist_ee_150_num, _hist_ee_150_den, bookScatter2D(2, 1, 3));
      divide(_hist_ee_240_num, _hist_ee_240_den, bookScatter2D(2, 1, 4));

      divide(_hist_ll_100_num, _hist_ll_100_den, bookScatter2D(3, 1, 1));
      divide(_hist_ll_125_num, _hist_ll_125_den, bookScatter2D(3, 1, 2));
      divide(_hist_ll_150_num, _hist_ll_150_den, bookScatter2D(3, 1, 3));
      divide(_hist_ll_240_num, _hist_ll_240_den, bookScatter2D(3, 1, 4));
    }

  private:
    /// Histograms
    Histo1D _hist_ee_100_num;
    Histo1D _hist_ee_125_num;
    Histo1D _hist_ee_150_num;
    Histo1D _hist_ee_240_num;

    Histo1D _hist_ee_100_den;
    Histo1D _hist_ee_125_den;
    Histo1D _hist_ee_150_den;
    Histo1D _hist_ee_240_den;

    Histo1D _hist_mm_100_num;
    Histo1D _hist_mm_125_num;
    Histo1D _hist_mm_150_num;
    Histo1D _hist_mm_240_num;

    Histo1D _hist_mm_100_den;
    Histo1D _hist_mm_125_den;
    Histo1D _hist_mm_150_den;
    Histo1D _hist_mm_240_den;

    Histo1D _hist_ll_100_num;
    Histo1D _hist_ll_125_num;
    Histo1D _hist_ll_150_num;
    Histo1D _hist_ll_240_num;

    Histo1D _hist_ll_100_den;
    Histo1D _hist_ll_125_den;
    Histo1D _hist_ll_150_den;
    Histo1D _hist_ll_240_den;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1122847);

}
