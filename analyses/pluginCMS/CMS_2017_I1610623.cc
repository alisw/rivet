#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"

#include "Rivet/AnalysisLoader.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Tools/RivetYODA.hh"

#include <iostream>

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2017_I1610623 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1610623);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs;
      WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 0*GeV, PID::MUON, 0*GeV, 1000000*GeV, 0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::TRACK, WFinder::TRANSMASS);
      //WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 0*GeV, PID::MUON, 0*GeV, 1000000*GeV, 0*GeV, 0.1, WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      addProjection(wfinder_mu, "WFinder_mu");

      // Define veto FS
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(wfinder_mu);
      vfs.addVetoPairId(PID::MUON);
      vfs.vetoNeutrinos();

      FastJets fastjets(vfs, FastJets::ANTIKT, 0.4);
      addProjection(fastjets, "Jets");

      //-------------
      _hist_Mult_exc      = bookHisto1D("d01-x01-y01");
      _hist_inc_WJetMult  = bookHisto1D("d02-x01-y01");

      //-------------
      _hist_JetPt1j = bookHisto1D("d03-x01-y01");
      _hist_JetPt2j = bookHisto1D("d04-x01-y01");
      _hist_JetPt3j = bookHisto1D("d05-x01-y01");
      _hist_JetPt4j = bookHisto1D("d06-x01-y01");

      //-------------
      _hist_JetRap1j = bookHisto1D("d07-x01-y01");
      _hist_JetRap2j = bookHisto1D("d08-x01-y01");
      _hist_JetRap3j = bookHisto1D("d09-x01-y01");
      _hist_JetRap4j = bookHisto1D("d10-x01-y01");

      //-------------
      _hist_Ht_1j = bookHisto1D("d11-x01-y01");
      _hist_Ht_2j = bookHisto1D("d12-x01-y01");
      _hist_Ht_3j = bookHisto1D("d13-x01-y01");
      _hist_Ht_4j = bookHisto1D("d14-x01-y01");

      //-------------
      _hist_dphij1mu_1j =bookHisto1D("d15-x01-y01");
      _hist_dphij2mu_2j =bookHisto1D("d16-x01-y01");
      _hist_dphij3mu_3j =bookHisto1D("d17-x01-y01");
      _hist_dphij4mu_4j =bookHisto1D("d18-x01-y01");

      //-------------
      _hist_dRmuj_1j =bookHisto1D("d19-x01-y01");

    }

    // define function used for filiing inc Njets histo
    void Fill(Histo1DPtr& _histJetMult, const double& weight, std::vector<FourMomentum>& finaljet_list){
      _histJetMult->fill(0, weight);
      for (size_t i=0 ; i<finaljet_list.size() ; ++i) {
        if (i==6) break;
        _histJetMult->fill(i+1, weight);  // inclusive multiplicity
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      const double weight = event.weight();
      const WFinder& wfinder_mu = applyProjection<WFinder>(event, "WFinder_mu");

      if (wfinder_mu.bosons().size() != 1) {
        vetoEvent;
      }

      if (wfinder_mu.bosons().size() == 1) {

        const FourMomentum lepton0 = wfinder_mu.constituentLepton().momentum();
        const FourMomentum neutrino = wfinder_mu.constituentNeutrino().momentum();
        double WmT = wfinder_mu.mT();

        if (WmT < 50.0*GeV) vetoEvent;

        double pt0 = lepton0.pT();
        double eta0 = lepton0.eta();

        if ( (fabs(eta0) > 2.4) || (pt0 < 25.0*GeV) ) vetoEvent;

        // Obtain the jets.
        vector<FourMomentum> finaljet_list;
        vector<FourMomentum> jet100_list;
        double HT = 0.0;

        // loop over jets in an event, pushback in finaljet_list collection
        foreach (const Jet& j, applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV)) {
          const double jrap = j.momentum().rap();
          const double jpt = j.momentum().pT();
          if ( (fabs(jrap) < 2.4) && (deltaR(lepton0, j.momentum()) > 0.4) ) {
            if(jpt > 30.0*GeV) {
              finaljet_list.push_back(j.momentum());
              HT += j.momentum().pT();
            }
            if(jpt > 100.0*GeV) {
              jet100_list.push_back(j.momentum());
            }
          }
        } // end looping over jets

        //---------------------- FILL HISTOGRAMS ------------------

        // Multiplicity exc plot.
        _hist_Mult_exc->fill(finaljet_list.size(), weight);

        // Multiplicity inc plot.
        Fill(_hist_inc_WJetMult, weight, finaljet_list);

        // dRmuj plot.
        double mindR(99999);
        if(jet100_list.size()>=1) {
          for (unsigned ji = 0; ji < jet100_list.size(); ji++){
            double dr_(9999);
            dr_ = fabs(deltaR(lepton0, jet100_list[ji]));
            if (dr_ < mindR){
              mindR = dr_;
            }
          }
          if(jet100_list[0].pT() > 300.0*GeV){
            _hist_dRmuj_1j->fill(mindR, weight);
          }
        }

        if(finaljet_list.size()>=1) {
          _hist_JetPt1j->fill(finaljet_list[0].pT(), weight);
          _hist_JetRap1j->fill(fabs(finaljet_list[0].rap()), weight);
          _hist_Ht_1j->fill(HT, weight);
          _hist_dphij1mu_1j->fill(deltaPhi(finaljet_list[0].phi(), lepton0.phi()), weight);
        }

        if(finaljet_list.size()>=2) {
          _hist_JetPt2j->fill(finaljet_list[1].pT(), weight);
          _hist_JetRap2j->fill(fabs(finaljet_list[1].rap()), weight);
          _hist_Ht_2j->fill(HT, weight);
          _hist_dphij2mu_2j->fill(deltaPhi(finaljet_list[1].phi(), lepton0.phi()), weight);
        }

        if(finaljet_list.size()>=3) {
          _hist_JetPt3j->fill(finaljet_list[2].pT(), weight);
          _hist_JetRap3j->fill(fabs(finaljet_list[2].rap()), weight);
          _hist_Ht_3j->fill(HT, weight);
          _hist_dphij3mu_3j->fill(deltaPhi(finaljet_list[2].phi(), lepton0.phi()), weight);
        }

        if(finaljet_list.size()>=4) {
          _hist_JetPt4j->fill(finaljet_list[3].pT(), weight);
          _hist_JetRap4j->fill(fabs(finaljet_list[3].rap()), weight);
          _hist_Ht_4j->fill(HT, weight);
          _hist_dphij4mu_4j->fill(deltaPhi(finaljet_list[3].phi(), lepton0.phi()), weight);
        }
      } // close the Wboson loop

    } //void loop


    /// Normalise histograms etc., after the run
    void finalize() {
      const double crossec = !std::isnan(crossSectionPerEvent()) ? crossSection() : 61526.7*picobarn;
      if (std::isnan(crossSectionPerEvent())){
        MSG_INFO("No valid cross-section given, using NNLO xsec calculated by FEWZ " << crossec/picobarn << " pb");
      }

      scale(_hist_Mult_exc, crossec/picobarn/sumOfWeights());
      scale(_hist_inc_WJetMult, crossec/picobarn/sumOfWeights());

      scale(_hist_JetPt1j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetPt2j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetPt3j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetPt4j, crossec/picobarn/sumOfWeights());

      scale(_hist_JetRap1j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetRap2j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetRap3j, crossec/picobarn/sumOfWeights());
      scale(_hist_JetRap4j, crossec/picobarn/sumOfWeights());

      scale(_hist_Ht_1j, crossec/picobarn/sumOfWeights());
      scale(_hist_Ht_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_Ht_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_Ht_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dphij1mu_1j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij2mu_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij3mu_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij4mu_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dRmuj_1j, crossec/picobarn/sumOfWeights());

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hist_Mult_exc;
    Histo1DPtr _hist_inc_WJetMult;

    Histo1DPtr _hist_JetPt1j;
    Histo1DPtr _hist_JetPt2j;
    Histo1DPtr _hist_JetPt3j;
    Histo1DPtr _hist_JetPt4j;

    Histo1DPtr _hist_JetRap1j;
    Histo1DPtr _hist_JetRap2j;
    Histo1DPtr _hist_JetRap3j;
    Histo1DPtr _hist_JetRap4j;

    Histo1DPtr _hist_Ht_1j;
    Histo1DPtr _hist_Ht_2j;
    Histo1DPtr _hist_Ht_3j;
    Histo1DPtr _hist_Ht_4j;

    Histo1DPtr _hist_dphij1mu_1j;
    Histo1DPtr _hist_dphij2mu_2j;
    Histo1DPtr _hist_dphij3mu_3j;
    Histo1DPtr _hist_dphij4mu_4j;

    Histo1DPtr _hist_dRmuj_1j;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1610623);


}
