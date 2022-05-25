// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Measurement of dimuon continuum in proton-nucleus collisions
  class E288_1981_I153009 : public Analysis {
  public:

    /// Constructor
    E288_1981_I153009()
      : Analysis("E288_1981_I153009")
    {   }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      Cut cut = Cuts::etaIn(-15.,15.);
      ZFinder zfinder(fs, cut, PID::MUON, 3.5*GeV, 30.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE );
      declare(zfinder, "ZFinder");

      // Book histograms
      // 400 GeV and y = 0.03
      Histo1DPtr dummy;
      _hist_pT_M_400.add(5., 6., book(dummy,9, 1, 1));
      _hist_pT_M_400.add(6., 7., book(dummy,9, 1, 2));
      _hist_pT_M_400.add(7., 8., book(dummy,9, 1, 3));
      _hist_pT_M_400.add(8., 9., book(dummy,9, 1, 4));
      _hist_pT_M_400.add(9., 10., book(dummy,9, 1, 5));
      _hist_pT_M_400.add(10.,11., book(dummy,9, 1, 6));
      _hist_pT_M_400.add(11.,12., book(dummy,9, 1, 7));
      _hist_pT_M_400.add(12.,13., book(dummy,9, 1, 8));
      _hist_pT_M_400.add(13.,14., book(dummy,9, 1, 9));

      int Nbin = 50;
      book(_h_m_DiMuon ,"DiMuon_mass",Nbin,0.0,30.0);
      book(_h_pT_DiMuon,"DiMuon_pT",Nbin,0.0,20.0);
      book(_h_y_DiMuon,"DiMuon_y",Nbin,-8.0, 8.0);
      book(_h_xF_DiMuon,"DiMuon_xF",Nbin, -1.5,  1.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double sqrts_tol = 10. ;
      if (!isCompatibleWithSqrtS(27.4, sqrts_tol)) {
        MSG_ERROR("Incorrect beam energy used: " << sqrtS()/GeV);
        throw Error("Unexpected sqrtS ! Only 27.4 GeV is supported");
      }

      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() >= 1) {

        double Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
        double Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
        double Zpl   = zfinder.bosons()[0].momentum().pz()/GeV;
        double Zy    = zfinder.bosons()[0].momentum().rapidity();
        //double ZE    = zfinder.bosons()[0].momentum().E();

        double xf = 2.*Zpl/sqrtS() ;
        _h_xF_DiMuon->fill(xf);
        _h_m_DiMuon->fill(Zmass/GeV);
        _h_pT_DiMuon->fill(Zpt);
        _h_y_DiMuon ->fill(Zy);
        double Zymin = -1.0;
        double Zymax = 1.03;
        double Z_y_width = Zymax - Zymin ;
        if ( Zy > Zymin && Zy < Zymax ) {
          // Edsigma^3/dp^3 = 2E/(pi*sqrts)dsigma/dx_F/dq_T^2 = 1/pi dsigma/dy/dq_T^2
          // normalisation of Zy bin width = Zwidth
          if (Zpt > 0) _hist_pT_M_400.fill(Zmass,Zpt, 1./2./Zpt/Z_y_width);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_DEBUG("Generator cross section [pb] = " << crossSection()/picobarn);
      _hist_pT_M_400.scale(crossSection()/femtobarn/(sumOfWeights() * M_PI), this);
    }

    //@}


    /// @name Histograms
    //@{
    BinnedHistogram _hist_pT_M_400,_hist_pT_M_300,_hist_pT_M_200;
    Histo1DPtr _h_m_DiMuon ;
    Histo1DPtr _h_pT_DiMuon;
    Histo1DPtr _h_y_DiMuon;
    Histo1DPtr _h_xF_DiMuon;
    Histo1DPtr _h_XXXX, _h_YYYY, _h_ZZZZ;
    Profile1DPtr _p_AAAA;
    CounterPtr _c_BBBB;
    //@}

  };


  RIVET_DECLARE_PLUGIN(E288_1981_I153009);

}
