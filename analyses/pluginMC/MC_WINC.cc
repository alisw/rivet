// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {



  /// @brief MC validation analysis for inclusive W events
  class MC_WINC : public Analysis {
  public:

    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_WINC);

    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
		  _dR=0.2;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
		  _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      FinalState fs;
      WFinder wfinder(fs, Cuts::abseta < 3.5 && Cuts::pT > 25*GeV, _lepton, 60.0*GeV, 100.0*GeV, 25.0*GeV, _dR);
      declare(wfinder, "WFinder");

      double sqrts = sqrtS()>0. ? sqrtS() : 14000.;
      book(_h_W_mass ,"W_mass", 50, 55.0, 105.0);
      book(_h_W_mT ,"W_mT", 40, 60.0, 100.0);
      book(_h_W_pT ,"W_pT", logspace(100, 1.0, 0.5*sqrts));
      book(_h_W_pT_peak ,"W_pT_peak", 25, 0.0, 125.0);
      book(_h_W_y ,"W_y", 40, -4.0, 4.0);
      book(_h_W_phi ,"W_phi", 25, 0.0, TWOPI);
      book(_h_Wplus_pT ,"Wplus_pT", logspace(25, 10.0, 0.5*sqrts));
      book(_h_Wminus_pT ,"Wminus_pT", logspace(25, 10.0, 0.5*sqrts));
      book(_h_lepton_pT ,"lepton_pT", logspace(100, 10.0, 0.25*sqrts));
      book(_h_lepton_eta ,"lepton_eta", 40, -4.0, 4.0);
      book(_htmp_dsigminus_deta ,"lepton_dsigminus_deta", 20, 0.0, 4.0);
      book(_htmp_dsigplus_deta  ,"lepton_dsigplus_deta", 20, 0.0, 4.0);

      book(_h_asym, "W_chargeasymm_eta");
      book(_h_asym_pT, "W_chargeasymm_pT");
    }



    /// Do the analysis
    void analyze(const Event & e) {

      const WFinder& wfinder = apply<WFinder>(e, "WFinder");
      if (wfinder.bosons().size() != 1) {
        vetoEvent;
      }

      double charge3_x_eta = 0;
      int charge3 = 0;
      FourMomentum emom;
      FourMomentum wmom(wfinder.bosons().front().momentum());
      _h_W_mass->fill(wmom.mass()/GeV);
      _h_W_mT->fill(wfinder.mT()/GeV);
      _h_W_pT->fill(wmom.pT()/GeV);
      _h_W_pT_peak->fill(wmom.pT()/GeV);
      _h_W_y->fill(wmom.rapidity());
      _h_W_phi->fill(wmom.phi());
      Particle l=wfinder.constituentLeptons()[0];
      _h_lepton_pT->fill(l.pT()/GeV);
      _h_lepton_eta->fill(l.eta());
      if (PID::charge3(l.pid()) != 0) {
        emom = l.momentum();
        charge3_x_eta = PID::charge3(l.pid()) * emom.eta();
        charge3 = PID::charge3(l.pid());
      }
      assert(charge3_x_eta != 0);
      assert(charge3!=0);
      if (emom.Et() > 30/GeV) {
        if (charge3_x_eta < 0) {
          _htmp_dsigminus_deta->fill(emom.eta());
        } else {
          _htmp_dsigplus_deta->fill(emom.eta());
        }
      }
      if (charge3 < 0) {
        _h_Wminus_pT->fill(wmom.pT()/GeV);
      } else {
        _h_Wplus_pT->fill(wmom.pT()/GeV);
      }
    }


    /// Finalize
    void finalize() {
      scale(_h_W_mass, crossSection()/sumOfWeights());
      scale(_h_W_mT, crossSection()/sumOfWeights());
      scale(_h_W_pT, crossSection()/sumOfWeights());
      scale(_h_W_pT_peak, crossSection()/sumOfWeights());
      scale(_h_W_y, crossSection()/sumOfWeights());
      scale(_h_W_phi, crossSection()/sumOfWeights());
      scale(_h_lepton_pT, crossSection()/sumOfWeights());
      scale(_h_lepton_eta, crossSection()/sumOfWeights());

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each Et region
      divide(*_htmp_dsigplus_deta - *_htmp_dsigminus_deta,
             *_htmp_dsigplus_deta + *_htmp_dsigminus_deta,
             _h_asym);

      // // W charge asymmetry vs. pTW: dsig+/dpT / dsig-/dpT
      divide(_h_Wplus_pT, _h_Wminus_pT,
             _h_asym_pT);

      scale(_h_Wplus_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_Wminus_pT, crossSection()/picobarn/sumOfWeights());
    }

    //@}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    //@{
    double _dR;
    PdgId _lepton;
    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_W_mass;
    Histo1DPtr _h_W_mT;
    Histo1DPtr _h_W_pT;
    Histo1DPtr _h_W_pT_peak;
    Histo1DPtr _h_W_y;
    Histo1DPtr _h_W_phi;
    Histo1DPtr _h_Wplus_pT;
    Histo1DPtr _h_Wminus_pT;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;

    Histo1DPtr _htmp_dsigminus_deta;
    Histo1DPtr _htmp_dsigplus_deta;

    Scatter2DPtr _h_asym;
    Scatter2DPtr _h_asym_pT;
    //@}

  };

  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WINC);
}
