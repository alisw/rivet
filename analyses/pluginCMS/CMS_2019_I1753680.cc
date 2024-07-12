// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Measurements of differential Z boson production cross sections in proton-proton collisions at 13 TeV
  class CMS_2019_I1753680 : public Analysis {
  public:
    
    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2019_I1753680);
    

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      // default to combined.
      _mode = 2;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;
      if ( getOption("LMODE") == "EMU" ) _mode = 2;

      // Initialise and register projections
      FinalState fs;
      Cut cut = Cuts::abseta < 2.4 && Cuts::pT > 25*GeV;

      ZFinder zeeFind(fs, cut, PID::ELECTRON, 76.1876*GeV, 106.1876*GeV, 0.1, ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES );
      declare(zeeFind, "ZeeFind");
      ZFinder zmmFind(fs, cut, PID::MUON    , 76.1876*GeV, 106.1876*GeV, 0.1, ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES );
      declare(zmmFind, "ZmmFind");
      
      // Book histograms
      book(_h_Zmm_absY          , 26, 1, 1);
      book(_h_Zee_absY          , 26, 1, 2);
      book(_h_Zll_absY          , 26, 1, 3);
      book(_h_Zmm_pt            , 27, 1, 1);
      book(_h_Zee_pt            , 27, 1, 2);
      book(_h_Zll_pt            , 27, 1, 3);
      book(_h_Zmm_phiStar       , 28, 1, 1);
      book(_h_Zee_phiStar       , 28, 1, 2);
      book(_h_Zll_phiStar       , 28, 1, 3);
      book(_h_Zll_pt_Y0         , 29, 1, 1);
      book(_h_Zll_pt_Y1         , 29, 1, 2);
      book(_h_Zll_pt_Y2         , 29, 1, 3);
      book(_h_Zll_pt_Y3         , 29, 1, 4);
      book(_h_Zll_pt_Y4         , 29, 1, 5);
      
      book(_h_Zll_pt_norm       , 30, 1, 1);
      book(_h_Zll_phiStar_norm  , 31, 1, 1);
      book(_h_Zll_absY_norm     , 32, 1, 1);
      book(_h_Zll_pt_Y0_norm    , 33, 1, 1);
      book(_h_Zll_pt_Y1_norm    , 33, 1, 2);
      book(_h_Zll_pt_Y2_norm    , 33, 1, 3);
      book(_h_Zll_pt_Y3_norm    , 33, 1, 4);
      book(_h_Zll_pt_Y4_norm    , 33, 1, 5);
      
    }
    

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const ZFinder& zeeFS = apply<ZFinder>(event, "ZeeFind");
      const ZFinder& zmumuFS = apply<ZFinder>(event, "ZmmFind");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      //event identification depending on mass window
      bool ee_event=false;
      bool mm_event=false;

      if (zees.size() == 1) { 
        ee_event = true; 
      }
      if (zmumus.size() == 1) { 
        mm_event = true; 
      }
      
      if (ee_event && _mode == 1)
        vetoEvent;
      if (mm_event && _mode == 0)
        vetoEvent;

      const Particles& theLeptons = ee_event ? zeeFS.constituents() : zmumuFS.constituents();
      const Particle& lminus = theLeptons[0].charge() < 0 ? theLeptons[0] : theLeptons[1];
      const Particle& lplus = theLeptons[0].charge() < 0 ? theLeptons[1] : theLeptons[0];

      //calculate phi*
      const double thetaStar = acos(tanh( 0.5 * (lminus.eta() - lplus.eta()) ));
      const double dPhi = M_PI - deltaPhi(lminus, lplus);
      const double phiStar = tan(0.5 * dPhi) * sin(thetaStar);
      
      const Particle& zcand = ee_event ? zees[0] : zmumus[0];

      if (ee_event) {
        _h_Zee_absY->fill(zcand.absrap());
        _h_Zee_pt->fill(zcand.pt());
        _h_Zee_phiStar->fill(phiStar);
      }
      else if (mm_event) {
        _h_Zmm_absY->fill(zcand.absrap());
        _h_Zmm_pt->fill(zcand.pt());
        _h_Zmm_phiStar->fill(phiStar);
      }

      _h_Zll_pt->fill(zcand.pt());
      _h_Zll_pt_norm->fill(zcand.pt());
      _h_Zll_phiStar->fill(phiStar);
      _h_Zll_phiStar_norm->fill(phiStar);
      _h_Zll_absY->fill(zcand.absrap());
      _h_Zll_absY_norm->fill(zcand.absrap());

      if      (zcand.absrap()<0.4) {
        _h_Zll_pt_Y0->fill(zcand.pt());
        _h_Zll_pt_Y0_norm->fill(zcand.pt());
      }
      else if (zcand.absrap()<0.8) {
        _h_Zll_pt_Y1->fill(zcand.pt());
        _h_Zll_pt_Y1_norm->fill(zcand.pt());
      }
      else if (zcand.absrap()<1.2) {
        _h_Zll_pt_Y2->fill(zcand.pt());
        _h_Zll_pt_Y2_norm->fill(zcand.pt());
      }
      else if (zcand.absrap()<1.6) {
        _h_Zll_pt_Y3->fill(zcand.pt());
        _h_Zll_pt_Y3_norm->fill(zcand.pt());
      }
      else if (zcand.absrap()<2.4) {
        _h_Zll_pt_Y4->fill(zcand.pt());
        _h_Zll_pt_Y4_norm->fill(zcand.pt());
      }

    }
    
    void normalizeToSum(Histo1DPtr hist) {
      double sum = 0.;
      for (size_t i = 0; i < hist->numBins(); ++i) {
        sum += hist->bin(i).height();
      }
      scale(hist, 1./sum);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/picobarn/sumOfWeights() : 1.0;
      
      scale(_h_Zmm_pt,      norm);
      scale(_h_Zmm_absY,    norm);
      scale(_h_Zmm_phiStar, norm);
      
      scale(_h_Zee_pt,      norm);
      scale(_h_Zee_absY,    norm);
      scale(_h_Zee_phiStar, norm);
      
      // when running in combined mode, need to average to get lepton xsec
      if (_mode == 2) norm /= 2.;
      
      scale(_h_Zll_pt,      norm);
      scale(_h_Zll_absY,    norm);
      scale(_h_Zll_phiStar, norm);
      scale(_h_Zll_pt_Y0,   norm);
      scale(_h_Zll_pt_Y1,   norm);
      scale(_h_Zll_pt_Y2,   norm);
      scale(_h_Zll_pt_Y3,   norm);
      scale(_h_Zll_pt_Y4,   norm);

      normalizeToSum(_h_Zll_pt_norm);
      normalizeToSum(_h_Zll_absY_norm);
      normalizeToSum(_h_Zll_phiStar_norm);
      normalizeToSum(_h_Zll_pt_Y0_norm);
      normalizeToSum(_h_Zll_pt_Y1_norm);
      normalizeToSum(_h_Zll_pt_Y2_norm);
      normalizeToSum(_h_Zll_pt_Y3_norm);
      normalizeToSum(_h_Zll_pt_Y4_norm);

    }

    //@}

  protected:

    size_t _mode;

    /// @name Histograms

  private:
    
    Histo1DPtr   _h_Zmm_pt, _h_Zmm_phiStar, _h_Zmm_absY;
    Histo1DPtr   _h_Zee_pt, _h_Zee_phiStar, _h_Zee_absY;

    Histo1DPtr   _h_Zll_pt, _h_Zll_phiStar, _h_Zll_absY;
    Histo1DPtr   _h_Zll_pt_Y0, _h_Zll_pt_Y1, _h_Zll_pt_Y2, _h_Zll_pt_Y3, _h_Zll_pt_Y4;

    Histo1DPtr   _h_Zll_pt_norm, _h_Zll_phiStar_norm, _h_Zll_absY_norm;
    Histo1DPtr   _h_Zll_pt_Y0_norm, _h_Zll_pt_Y1_norm, _h_Zll_pt_Y2_norm, _h_Zll_pt_Y3_norm, _h_Zll_pt_Y4_norm;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2019_I1753680);


}
