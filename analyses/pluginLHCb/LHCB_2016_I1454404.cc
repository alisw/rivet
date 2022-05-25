// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Measurement of forward W and Z boson production with jets in pp collisions at 8 TeV
  class LHCB_2016_I1454404 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2016_I1454404);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      _mode = 0;
      string mode = getOption("MODE");
      if (mode == "ALL" ) _mode = 0;
      else if (mode == "WpJET") _mode = 1;
      else if (mode == "WmJET") _mode = 2;
      else if (mode == "ZJET") _mode = 3;
      else if (mode == "WJET") _mode = 4;


      const Cut muSel = Cuts::eta >= 2.0 && Cuts::eta <= 4.5 && Cuts::pT > 20*GeV;

      // Z boson
      ZFinder zfinder(FinalState(), muSel, PID::MUON, 60*GeV, 120*GeV);
      declare(zfinder, "ZFinder");

      // W boson
      WFinder wfinder(FinalState(), muSel, PID::MUON, 0*GeV ,500*GeV ,0.);
      declare(wfinder, "WFinder");

      // Jet Z
      FastJets jetproZ(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.5);
      declare(jetproZ, "JetsZ");

      // Jet W
      FastJets jetproW(wfinder.remainingFinalState(), FastJets::ANTIKT, 0.5);
      declare(jetproW, "JetsW");

      // Book histograms
      /////////
      if (_mode == 0 || _mode == 1 || _mode == 4) {
        book(_h_wpj, 1, 1, 1);
        book(_h_eta_wpj, 4, 1, 1);
        book(_h_etaj_wpj, 5, 1, 1);
        book(_h_ptj_wpj, 6, 1, 1);
      }
      /////////
      if (_mode == 0 || _mode == 2 || _mode == 4) {
        book(_h_wmj, 1, 1, 2);
        book(_h_eta_wmj, 4, 1, 2);
        book(_h_etaj_wmj, 5, 1, 2);
        book(_h_ptj_wmj, 6, 1, 2);
      }

      /////////
      if (_mode == 0 || _mode == 3) {
        book(_h_zj, 1, 1, 3);
        book(_h_yz_zj, 7, 1, 1);
        book(_h_etaj_zj, 8, 1, 1);
        book(_h_ptj_zj, 9, 1, 1);
        book(_h_dphi_zj, 10, 1, 1);
      }

      if (_mode == 0 ){
        book(_h_rwz, 2,1,1);
        book(_h_rwpz, 2,1,2);
        book(_h_rwmz, 2,1,3);
      }

      if (_mode == 0 || _mode == 4){
        book(_h_rwpm, 2,1,4);
        book(_h_aw, 3,1,1);
        // this is a temporary histogram to construct rwz later
        book(_h_wj, "_temp_wj", refData(1,1,1));
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Cut jetSel = Cuts::eta >= 2.2 && Cuts::eta <= 4.2 && Cuts::pT > 20*GeV;
      if (_mode == 0 || _mode == 3) {

        //////////////////////////////////////////////////////////
        ///////////////ZFinder Muon //////////////////////////////
        //////////////////////////////////////////////////////////

        const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
        if (zfinder.bosons().size() ==1){
          const Particles muon = zfinder.constituentLeptons(); //zfinder.constituents()?
          const Particles Z = zfinder.bosons();
          const FourMomentum Zmom = Z[0].momentum();
          const Jets jetsZ = apply<FastJets>(event, "JetsZ").jetsByPt(jetSel);
          const Jets cleanedJetsZ = filter_discard(jetsZ, [&](const Jet& j) {return any(muon, deltaRLess(j, 0.5)); });

          if (cleanedJetsZ.size() > 0 && cleanedJetsZ.at(0).pT() > 20*GeV) {
            const double yZ = Zmom.rap(); //histogram 7
            const double etaj = cleanedJetsZ[0].eta(); //histogram 8
            const double ptj  = cleanedJetsZ[0].pT()/GeV; //histogram 9
            double dphi_tmp = abs(Zmom.phi() - cleanedJetsZ[0].phi());
            const double dphi = dphi_tmp < Rivet::pi ? dphi_tmp : Rivet::twopi - dphi_tmp;
            _h_zj->fill(sqrtS()/GeV);
            _h_dphi_zj->fill(dphi);
            _h_yz_zj->fill(yZ); // boson rapidity vs diff cross section
            _h_etaj_zj->fill(etaj); // jet pseudorapidity vs diff cross section
            _h_ptj_zj->fill(ptj); //jet transverse momentum vs diff cross section
          }
        }
      }

      if (_mode == 0 || _mode == 1 || _mode == 2  || _mode == 4) {
        //////////////////////////////////////////////////////////
        ///////////////WFinder Muon //////////////////////////////
        //////////////////////////////////////////////////////////

        const WFinder& wfinder = apply<WFinder>(event, "WFinder");
        if (wfinder.bosons().size() == 1) {
          const Particles Muons = wfinder.constituentLeptons();
          const FourMomentum muonmom = Muons[0].momentum();
          const Jets jetsW = apply<FastJets>(event, "JetsW").jetsByPt(jetSel);

          const Jets cleanedJetsW = filter_discard(jetsW, [&](const Jet& j) {return any(Muons, deltaRLess(j, 0.5)); });

          if (cleanedJetsW.size() > 0 && cleanedJetsW.at(0).pT() > 20*GeV) {
            const double etaj = cleanedJetsW[0].eta(); //histogram 5
            const double etamu = muonmom.eta(); //histogram 4
            if( (_mode == 0 || _mode == 1 || _mode == 4) && Muons[0].charge() > 0) {
              //fill with W related analysis
              if (_mode != 1 ) _h_wj->fill(sqrtS()/GeV); // don't need this for single charge case
              _h_wpj->fill(sqrtS()/GeV);
              _h_eta_wpj->fill(etamu); // W+ Jet  muon pseudorapidity vs diff cross section
              _h_etaj_wpj->fill(etaj); // W+ Jet jet pseudorapidity vs diff cross section
              _h_ptj_wpj->fill(cleanedJetsW[0].pT()/GeV); // W+ Jet jet transverse momentum vs diff cross section
            }
            else if( (_mode == 0 || _mode == 2 || _mode == 4 ) && Muons[0].charge() < 0) {
              //fill with W related analysis
              if (_mode != 2) _h_wj->fill(sqrtS()/GeV);
              _h_wmj->fill(sqrtS()/GeV); // don't need this for single charge case
              _h_eta_wmj->fill(etamu); // W- Jet  muon pseudorapidity vs diff cross section
              _h_etaj_wmj->fill(etaj); // W- Jet jet pseudorapidity vs diff cross section
              _h_ptj_wmj->fill(cleanedJetsW[0].pT()/GeV); // W+ Jet jet transverse momentum vs diff cross section
            }
          }
        }

      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      double scalefactor = crossSection()/picobarn/sumOfWeights();

      if(_mode == 0 || _mode == 1 || _mode == 4) {
        scale({_h_wpj, _h_eta_wpj, _h_etaj_wpj, _h_ptj_wpj}, scalefactor);
      }
      if(_mode == 0 || _mode == 2 || _mode == 4) {
        scale({_h_wmj, _h_eta_wmj, _h_etaj_wmj, _h_ptj_wmj}, scalefactor);
      }
      if(_mode == 0 || _mode == 3) {
        scale({_h_zj, _h_yz_zj, _h_etaj_zj, _h_ptj_zj, _h_dphi_zj}, scalefactor);
      }
      if (_mode == 0 ) {
        scale(_h_wj, scalefactor); // need to scale this for consistency
        divide(_h_wpj, _h_zj, _h_rwpz);
        divide(_h_wmj, _h_zj, _h_rwmz);
        divide(_h_wj, _h_zj, _h_rwz);

      }
      if (_mode == 0 || _mode == 4) {
        divide(_h_wpj, _h_wmj, _h_rwpm);
        asymm(_h_wpj, _h_wmj, _h_aw);
      }

    }

  protected:

    Log& log = getLog(); // in Analysis or Projection
    size_t _mode;

    /// histograms
    Histo1DPtr _h_wpj, _h_wmj, _h_wj, _h_zj;
    Scatter2DPtr _h_rwz, _h_rwpz, _h_rwmz, _h_rwpm, _h_aw;
    Histo1DPtr _h_eta_wpj, _h_eta_wmj, _h_etaj_wpj, _h_etaj_wmj, _h_ptj_wpj, _h_ptj_wmj;
    Histo1DPtr _h_yz_zj, _h_etaj_zj, _h_ptj_zj, _h_dphi_zj;
  };


  RIVET_DECLARE_PLUGIN(LHCB_2016_I1454404);

}
