// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for W^+[enu]W^-[munu] events
  class MC_WWINC : public Analysis {
  public:

    /// Default constructor
    MC_WWINC()
      : Analysis("MC_WWINC")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      WFinder wenufinder(fs, Cuts::abseta < 3.5 && Cuts::pT > 25*GeV, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wenufinder, "WenuFinder");

      VetoedFinalState wmnuinput;
      wmnuinput.addVetoOnThisFinalState(wenufinder);
      WFinder wmnufinder(wmnuinput, Cuts::abseta < 3.5 && Cuts::pT > 25*GeV, PID::MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wmnufinder, "WmnuFinder");

      // properties of the pair momentum
      double sqrts = sqrtS()>0. ? sqrtS() : 14000.;
      book(_h_WW_pT ,"WW_pT", logspace(100, 1.0, max(1.1,0.5*sqrts)));
      book(_h_WW_pT_peak ,"WW_pT_peak", 25, 0.0, 25.0);
      book(_h_WW_eta ,"WW_eta", 40, -7.0, 7.0);
      book(_h_WW_phi ,"WW_phi", 25, 0.0, TWOPI);
      book(_h_WW_m ,"WW_m", logspace(100, 150.0, 180.0+0.25*sqrts));

      // correlations between the WW
      book(_h_WW_dphi ,"WW_dphi", 25, 0.0, PI);  /// @todo non-linear?
      book(_h_WW_deta ,"WW_deta", 25, -7.0, 7.0);
      book(_h_WW_dR ,"WW_dR", 25, 0.5, 7.0);
      book(_h_WW_dpT ,"WW_dpT", logspace(100, 1.0, max(1.1,0.5*sqrts)));
      book(_h_WW_costheta_planes ,"WW_costheta_planes", 25, -1.0, 1.0);

      /// @todo fuer WW: missing ET

      // properties of the W bosons
      book(_h_W_pT ,"W_pT", logspace(100, 10.0, max(11.,0.25*sqrts)));
      book(_h_W_eta ,"W_eta", 70, -7.0, 7.0);

      // properties of the leptons
      book(_h_Wl_pT ,"Wl_pT", logspace(100, 30.0, max(31., 0.1*sqrts)));
      book(_h_Wl_eta ,"Wl_eta", 40, -3.5, 3.5);

      // correlations between the opposite charge leptons
      book(_h_WeWm_dphi ,"WeWm_dphi", 25, 0.0, PI);
      book(_h_WeWm_deta ,"WeWm_deta", 25, -5.0, 5.0);
      book(_h_WeWm_dR ,"WeWm_dR", 25, 0.5, 5.0);
      book(_h_WeWm_m ,"WeWm_m", 100, 0.0, 300.0);
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = 1.0;

      const WFinder& wenufinder = apply<WFinder>(e, "WenuFinder");
      if (wenufinder.bosons().size()!=1) {
        vetoEvent;
      }

      const WFinder& wmnufinder = apply<WFinder>(e, "WmnuFinder");
      if (wmnufinder.bosons().size()!=1) {
        vetoEvent;
      }

      FourMomentum wenu(wenufinder.bosons()[0].momentum());
      FourMomentum wmnu(wmnufinder.bosons()[0].momentum());
      FourMomentum ww(wenu+wmnu);
      // find leptons
      FourMomentum ep=wenufinder.constituentLeptons()[0].momentum();
      FourMomentum enu=wenufinder.constituentNeutrinos()[0].momentum();
      FourMomentum mm=wmnufinder.constituentLeptons()[0].momentum();
      FourMomentum mnu=wmnufinder.constituentNeutrinos()[0].momentum();

      _h_WW_pT->fill(ww.pT(),weight);
      _h_WW_pT_peak->fill(ww.pT(),weight);
      _h_WW_eta->fill(ww.eta(),weight);
      _h_WW_phi->fill(ww.phi(),weight);
      double mww2=ww.mass2();
      if (mww2>0.0) _h_WW_m->fill(sqrt(mww2), weight);

      _h_WW_dphi->fill(mapAngle0ToPi(wenu.phi()-wmnu.phi()), weight);
      _h_WW_deta->fill(wenu.eta()-wmnu.eta(), weight);
      _h_WW_dR->fill(deltaR(wenu,wmnu), weight);
      _h_WW_dpT->fill(fabs(wenu.pT()-wmnu.pT()), weight);

      Vector3 crossWenu = ep.p3().cross(enu.p3());
      Vector3 crossWmnu = mm.p3().cross(mnu.p3());
      double costheta = crossWenu.dot(crossWmnu)/crossWenu.mod()/crossWmnu.mod();
      _h_WW_costheta_planes->fill(costheta, weight);

      _h_W_pT->fill(wenu.pT(),weight);
      _h_W_pT->fill(wmnu.pT(),weight);
      _h_W_eta->fill(wenu.eta(),weight);
      _h_W_eta->fill(wmnu.eta(),weight);

      _h_Wl_pT->fill(ep.pT(), weight);
      _h_Wl_pT->fill(mm.pT(), weight);
      _h_Wl_eta->fill(ep.eta(), weight);
      _h_Wl_eta->fill(mm.eta(), weight);

      _h_WeWm_dphi->fill(mapAngle0ToPi(ep.phi()-mm.phi()), weight);
      _h_WeWm_deta->fill(ep.eta()-mm.eta(), weight);
      _h_WeWm_dR->fill(deltaR(ep,mm), weight);
      double m2=FourMomentum(ep+mm).mass2();
      if (m2 < 0) m2 = 0.0;
      _h_WeWm_m->fill(sqrt(m2), weight);
    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/picobarn/sumOfWeights();
      scale(_h_WW_pT, norm);
      scale(_h_WW_pT_peak, norm);
      scale(_h_WW_eta, norm);
      scale(_h_WW_phi, norm);
      scale(_h_WW_m, norm);
      scale(_h_WW_dphi, norm);
      scale(_h_WW_deta, norm);
      scale(_h_WW_dR, norm);
      scale(_h_WW_dpT, norm);
      scale(_h_WW_costheta_planes, norm);
      scale(_h_W_pT, norm);
      scale(_h_W_eta, norm);
      scale(_h_Wl_pT, norm);
      scale(_h_Wl_eta, norm);
      scale(_h_WeWm_dphi, norm);
      scale(_h_WeWm_deta, norm);
      scale(_h_WeWm_dR, norm);
      scale(_h_WeWm_m, norm);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_WW_pT;
    Histo1DPtr _h_WW_pT_peak;
    Histo1DPtr _h_WW_eta;
    Histo1DPtr _h_WW_phi;
    Histo1DPtr _h_WW_m;
    Histo1DPtr _h_WW_dphi;
    Histo1DPtr _h_WW_deta;
    Histo1DPtr _h_WW_dR;
    Histo1DPtr _h_WW_dpT;
    Histo1DPtr _h_WW_costheta_planes;
    Histo1DPtr _h_W_pT;
    Histo1DPtr _h_W_eta;
    Histo1DPtr _h_Wl_pT;
    Histo1DPtr _h_Wl_eta;
    Histo1DPtr _h_WeWm_dphi;
    Histo1DPtr _h_WeWm_deta;
    Histo1DPtr _h_WeWm_dR;
    Histo1DPtr _h_WeWm_m;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WWINC);

}
