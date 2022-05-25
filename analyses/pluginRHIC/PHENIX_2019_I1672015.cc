// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Drell Yan production at low masses at $\sqrt{s} = 200$ GeV
  class PHENIX_2019_I1672015 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PHENIX_2019_I1672015);

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      declare(fs, "FS");
      Cut cut = Cuts::etaIn(-10.,10.);
      ZFinder zfinder(fs, cut, PID::MUON, 4.0*GeV, 100.0*GeV, 0.1, ZFinder::ClusterPhotons::NONE );
      declare(zfinder, "ZFinder");

      // Book histograms
      book(_h_pT ,1, 1, 1);
      book(_h_mass,2, 1, 1);
      int Nbin = 50;
      book(_h_m_DiMuon,"DiMuon_mass",Nbin,0.0,30.0);
      book(_h_pT_DiMuon,"DiMuon_pT",Nbin,0.0,20.0);
      book(_h_y_DiMuon,"DiMuon_y",Nbin,-8.0, 8.0);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double sqrts_tol = 10. ;
      if (!isCompatibleWithSqrtS(200., sqrts_tol)) {
        MSG_ERROR("Incorrect beam energy used: " << sqrtS()/GeV);
        throw Error("Unexpected sqrtS ! Only 200 GeV is supported");
      }
      const ZFinder& zfinder = applyProjection<ZFinder>(event, "ZFinder");
      if (zfinder.particles().size() < 1) vetoEvent;
      double mass = zfinder.bosons()[0].momentum().mass()/GeV;
      double pt_DY   = zfinder.bosons()[0].momentum().pT()/GeV;
      double y_DY    = abs(zfinder.bosons()[0].momentum().rapidity());

      _h_m_DiMuon->fill(mass/GeV);
      _h_pT_DiMuon->fill(pt_DY);
      _h_y_DiMuon ->fill(y_DY);
      if (inRange( y_DY, 1.2, 2.2) && inRange(mass, 4.8, 8.2) && pt_DY > 0)
        _h_pT->fill(pt_DY, 1./2./pt_DY/M_PI);
      if ( y_DY >= 1.2 && y_DY <= 2.2)
        _h_mass->fill(mass);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Divide by 2 because of rapidity range
      normalize(_h_m_DiMuon);
      normalize(_h_pT_DiMuon);
      normalize(_h_y_DiMuon);
      scale(_h_pT, crossSection()/(sumOfWeights())/2.);
      scale(_h_mass, crossSection()/(sumOfWeights())/2.);

    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pT, _h_y, _h_mass, _h_ZZZZ;
    Histo1DPtr _h_m_DiMuon;
    Histo1DPtr _h_pT_DiMuon;
    Histo1DPtr _h_y_DiMuon;
    Histo1DPtr _h_xF_DiMuon;
    ///@}

  };


  RIVET_DECLARE_PLUGIN(PHENIX_2019_I1672015);

}
