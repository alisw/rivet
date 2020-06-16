// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Dimuon production in proton-copper collisions at 38.8 GeV
  class E605_1991_I302822 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(E605_1991_I302822);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      Cut cut = Cuts::etaIn(-10.,10.);
      ZFinder zfinder(fs, cut, PID::MUON, 4.0*GeV, 100.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE );
      declare(zfinder, "ZFinder");

      // Book histograms in mass ranges (measurement is not normalised to mass range)
      Histo1DPtr dummy;
      book(_hist_pT_M_78,17, 1, 1);
      book(_hist_pT_M_89,18, 1, 1);
      book(_hist_pT_M_1011,19, 1, 1);
      book(_hist_pT_M_1113,20, 1, 1);
      book(_hist_pT_M_1318,21, 1, 1);

      int Nbin = 50;
      book(_h_m_DiMuon,"DiMuon_mass",Nbin,0.0,30.0);
      book(_h_pT_DiMuon,"DiMuon_pT",Nbin,0.0,20.0);
      book(_h_y_DiMuon,"DiMuon_y",Nbin,-8.0, 8.0);
      book(_h_xF_DiMuon,"DiMuon_xF",Nbin, -1.5,  1.5);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double sqrts_tol = 10. ;
      if (!fuzzyEquals(sqrtS()/GeV, 38.8, sqrts_tol)) {
        MSG_ERROR("Incorrect beam energy used: " << sqrtS()/GeV);
        throw Error("Unexpected sqrtS ! Only 38.8 GeV is supported");
      }

      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() >= 1) {

        double Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
        double Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
        double Zpl   = zfinder.bosons()[0].momentum().pz()/GeV;
        double Zy    = zfinder.bosons()[0].momentum().rapidity();
        double ZE    = zfinder.bosons()[0].momentum().E();

        double xf = 2.*Zpl/sqrtS();
        _h_xF_DiMuon->fill(xf);
        _h_m_DiMuon->fill(Zmass/GeV);
        _h_pT_DiMuon->fill(Zpt);
        _h_y_DiMuon ->fill(Zy);

        if (xf > -0.1 && xf < 0.2 && Zpt > 0) {
          if (Zmass > 7. && Zmass < 8.) _hist_pT_M_78->fill(Zpt, 1./2./Zpt *2.*ZE/sqrtS());
          if (Zmass > 8. && Zmass < 9.) _hist_pT_M_89->fill(Zpt, 1./2./Zpt *2.*ZE/sqrtS());
          if (Zmass > 10.5 && Zmass < 11.5) _hist_pT_M_1011->fill(Zpt, 1./2./Zpt *2.*ZE/sqrtS());
          if (Zmass > 11.5 && Zmass < 13.5) _hist_pT_M_1113->fill(Zpt, 1./2./Zpt *2.*ZE/sqrtS());
          if (Zmass > 13.5 && Zmass < 18.)  _hist_pT_M_1318->fill(Zpt, 1./2./Zpt *2.*ZE/sqrtS());
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_DEBUG(" Generator cross section [pb] " << crossSection()/picobarn);

      normalize(_h_m_DiMuon);
      normalize(_h_pT_DiMuon);
      normalize(_h_xF_DiMuon);
      normalize(_h_y_DiMuon);

      // normalisation of xF bin width = 0.3
      const double scalefactor=crossSection()/picobarn/(sumOfWeights() * M_PI *0.3 );
      scale(_hist_pT_M_78,scalefactor);
      scale(_hist_pT_M_89,scalefactor);
      scale(_hist_pT_M_1011,scalefactor);
      scale(_hist_pT_M_1113,scalefactor);
      scale(_hist_pT_M_1318,scalefactor);
    }

    //@}


    /// @name Histograms
    //@{
    BinnedHistogram _hist_pT_M;
    Histo1DPtr _h_m_DiMuon ;
    Histo1DPtr _h_pT_DiMuon;
    Histo1DPtr _h_y_DiMuon;
    Histo1DPtr _h_xF_DiMuon;
    Histo1DPtr _hist_pT_M_78,_hist_pT_M_89,_hist_pT_M_1011,_hist_pT_M_1113,_hist_pT_M_1318;
    //@}

  };


  DECLARE_RIVET_PLUGIN(E605_1991_I302822);

}
