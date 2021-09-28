// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Drell Yan measurements $pp \to \mu^+\mu^- +X $ at $\sqrt{s} = 44$ and $62$ GeV at CERN ISR
  class R209_1982_I168182 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(R209_1982_I168182);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs;
      declare(fs, "FS");
      Cut cut = Cuts::etaIn(-10.,10.);
      ZFinder zfinder(fs, cut, PID::MUON, 3.5*GeV, 18.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE );
      declare(zfinder, "ZFinder");

      // Book histograms
      if (fuzzyEquals(sqrtS()/GeV, 62., sqrts_tol)) {
        MSG_DEBUG("R209: running with 62: " << sqrtS()/GeV);
        book(_hist_M,1, 1, 1);
        book(_hist_pT ,2, 1, 1);
      }
      else if (fuzzyEquals(sqrtS()/GeV, 44., sqrts_tol)) {
        MSG_DEBUG("R209: running with 44: " << sqrtS()/GeV);
        book(_hist_M,1, 1, 2);
      }
      int Nbin = 50;
      book(_h_m_DiMuon,"DiMuon_mass",Nbin,0.0,30.0);
      book(_h_pT_DiMuon,"DiMuon_pT",Nbin,0.0,20.0);
      book(_h_y_DiMuon,"DiMuon_y",Nbin,-8.0, 8.0);
      book(_h_xF_DiMuon,"DiMuon_xF",Nbin, -1.5,  1.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() >= 1) {

        double Zmass = zfinder.bosons()[0].momentum().mass()/GeV;
        double Zpt   = zfinder.bosons()[0].momentum().pT()/GeV;
        double Zpl   = zfinder.bosons()[0].momentum().pz()/GeV;
        double Zy    = zfinder.bosons()[0].momentum().rapidity();
        double xf = 2.*Zpl/sqrtS() ;

        _h_xF_DiMuon->fill(xf);
        _h_m_DiMuon->fill(Zmass/GeV);
        _h_pT_DiMuon->fill(Zpt);
        _h_y_DiMuon->fill(Zy);
        if (fuzzyEquals(sqrtS()/GeV, 62, sqrts_tol)) {
          if (Zmass > 0) _hist_M->fill(Zmass);
          if (Zmass > 5. && Zmass < 8.) {
            if (Zpt > 0) _hist_pT->fill(Zpt,1./2./Zpt);
          }
        }
        else if (fuzzyEquals(sqrtS()/GeV, 44, sqrts_tol)) {
          if (Zmass > 0) _hist_M->fill(Zmass);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_m_DiMuon);
      normalize(_h_pT_DiMuon);
      normalize(_h_xF_DiMuon);
      normalize(_h_y_DiMuon);
      scale(_hist_pT,crossSection()/nanobarn/(sumOfWeights()));
      scale(_hist_M,crossSection()/nanobarn/(sumOfWeights()));
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _hist_pT, _hist_M ;
    Histo1DPtr _h_m_DiMuon ;
    Histo1DPtr _h_pT_DiMuon;
    Histo1DPtr _h_y_DiMuon;
    Histo1DPtr _h_xF_DiMuon;
    ///@}

    /// Energy comparison tolerance
    const double sqrts_tol = 0.1;

  };


  DECLARE_RIVET_PLUGIN(R209_1982_I168182);

}
