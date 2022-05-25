// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class D0_2010_S8821313 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2010_S8821313);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections
      FinalState fs;
      Cut cuts = (Cuts::abseta < 1.1 || Cuts::absetaIn( 1.5,  3.0)) && Cuts::pT > 20*GeV;
      ZFinder zfinder_ee(fs, cuts, PID::ELECTRON, 70*GeV, 110*GeV, 0.2, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
      declare(zfinder_ee, "zfinder_ee");
      ZFinder zfinder_mm(fs, Cuts::abseta < 2 && Cuts::pT > 15*GeV, PID::MUON, 70*GeV, 110*GeV, 0.0, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::NO);
      declare(zfinder_mm, "zfinder_mm");

      /// Book histograms here
      {Histo1DPtr tmp; _h_phistar_ee.add(0.0, 1.0, book(tmp, 1, 1, 1));}
      {Histo1DPtr tmp; _h_phistar_ee.add(1.0, 2.0, book(tmp, 1, 1, 2));}
      {Histo1DPtr tmp; _h_phistar_ee.add(2.0, 10.0,book(tmp, 1, 1, 3));}
      {Histo1DPtr tmp; _h_phistar_mm.add(0.0, 1.0, book(tmp, 2, 1, 1));}
      {Histo1DPtr tmp; _h_phistar_mm.add(1.0, 2.0, book(tmp, 2, 1, 2));}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      const ZFinder& zfinder_ee = apply<ZFinder>(event, "zfinder_ee");
      if (zfinder_ee.bosons().size() == 1) {
        Particles ee = zfinder_ee.constituents();
        std::sort(ee.begin(), ee.end(), cmpMomByPt);
        const FourMomentum& eminus = PID::charge3(ee[0].pid()) < 0 ? ee[0].momentum() : ee[1].momentum();
        const FourMomentum& eplus  = PID::charge3(ee[0].pid()) < 0 ? ee[1].momentum() : ee[0].momentum();
        double phi_acop = M_PI - mapAngle0ToPi(eminus.phi() - eplus.phi());
        double costhetastar = tanh((eminus.eta() - eplus.eta())/2);
        double sin2thetastar = 1 - sqr(costhetastar);
        if (sin2thetastar < 0) sin2thetastar = 0;
        const double phistar = tan(phi_acop/2) * sqrt(sin2thetastar);
        const FourMomentum& zmom = zfinder_ee.bosons()[0].momentum();
        _h_phistar_ee.fill(zmom.rapidity(), phistar, weight);
      }

      const ZFinder& zfinder_mm = apply<ZFinder>(event, "zfinder_mm");
      if (zfinder_mm.bosons().size() == 1) {
        Particles mm = zfinder_mm.constituents();
        std::sort(mm.begin(), mm.end(), cmpMomByPt);
        const FourMomentum& mminus = PID::charge3(mm[0].pid()) < 0 ? mm[0].momentum() : mm[1].momentum();
        const FourMomentum& mplus  = PID::charge3(mm[0].pid()) < 0 ? mm[1].momentum() : mm[0].momentum();
        double phi_acop = M_PI - mapAngle0ToPi(mminus.phi() - mplus.phi());
        double costhetastar = tanh((mminus.eta() - mplus.eta())/2);
        double sin2thetastar = 1 - sqr(costhetastar);
        if (sin2thetastar < 0) sin2thetastar = 0;
        const double phistar = tan(phi_acop/2) * sqrt(sin2thetastar);
        const FourMomentum& zmom = zfinder_mm.bosons()[0].momentum();
        _h_phistar_mm.fill(zmom.rapidity(), phistar, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (Histo1DPtr hist : _h_phistar_ee.histos()) normalize(hist);
      for (Histo1DPtr hist : _h_phistar_mm.histos()) normalize(hist);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    BinnedHistogram _h_phistar_ee;
    BinnedHistogram _h_phistar_mm;
    //@}


  };



  RIVET_DECLARE_ALIASED_PLUGIN(D0_2010_S8821313, D0_2010_I871787);

}
