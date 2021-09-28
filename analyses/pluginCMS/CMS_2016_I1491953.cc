#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  /// @brief Differential cross sections for associated production of a W boson and jets at 8 TeV
  class CMS_2016_I1491953 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1491953);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState fs;
      WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 0*GeV, PID::MUON, 0*GeV, 1000000*GeV,
                         0*GeV, 0.1, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::YES, WFinder::MassWindow::MT);
      declare(wfinder_mu, "WFinder_mu");

      // Define veto FS
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(wfinder_mu);
      vfs.addVetoPairId(PID::MUON);
      vfs.vetoNeutrinos();

      FastJets fastjets(vfs, FastJets::ANTIKT, 0.5);
      declare(fastjets, "Jets");

      book(_hist_Mult_exc      ,1,1,1);
      book(_hist_inc_WJetMult  ,2,1,1);

      book(_hist_addJetPt1j ,3,1,1);
      book(_hist_addJetPt2j ,4,1,1);
      book(_hist_addJetPt3j ,5,1,1);
      book(_hist_addJetPt4j ,6,1,1);

      book(_hist_addHt_1j ,7,1,1);
      book(_hist_addHt_2j ,8,1,1);
      book(_hist_addHt_3j ,9,1,1);
      book(_hist_addHt_4j ,10,1,1);

      book(_hist_diJetPt_2j ,11,1,1);
      book(_hist_diJetPt_3j ,12,1,1);
      book(_hist_diJetPt_4j ,13,1,1);

      book(_hist_dijetM_2j ,14,1,1);
      book(_hist_dijetM_3j ,15,1,1);
      book(_hist_dijetM_4j ,16,1,1);

      book(_hist_Jeteta1j ,17,1,1);
      book(_hist_Jeteta2j ,18,1,1);
      book(_hist_Jeteta3j ,19,1,1);
      book(_hist_Jeteta4j ,20,1,1);

      book(_hist_dyj1j2_2j ,21,1,1);
      book(_hist_dyj1j2_3j ,22,1,1);
      book(_hist_dyj1j2_4j ,23,1,1);

      book(_hist_dyj1j3_3j ,24,1,1);
      book(_hist_dyj2j3_3j ,25,1,1);

      book(_hist_dyjFjB_2j ,26,1,1);
      book(_hist_dyjFjB_3j ,27,1,1);
      book(_hist_dyjFjB_4j ,28,1,1);

      book(_hist_dphij1j2_2j ,29,1,1);
      book(_hist_dphijFjB_2j ,30,1,1);
      book(_hist_dRj1j2_2j ,31,1,1);

      book(_hist_dphij1mu_1j ,32,1,1);
      book(_hist_dphij2mu_2j ,33,1,1);
      book(_hist_dphij3mu_3j ,34,1,1);
      book(_hist_dphij4mu_4j ,35,1,1);

      book(_hist_MeanNJht_1j     ,36,1,1);
      book(_hist_MeanNJht_2j     ,37,1,1);
      book(_hist_MeanNJdyj1j2_2j ,38,1,1);
      book(_hist_MeanNJdyjFjB_2j ,39,1,1);

    }


    // Define function used for filiing inc Njets histo
    void _fill(Histo1DPtr& _histJetMult, vector<FourMomentum>& finaljet_list) {
      _histJetMult->fill(0);
      for (size_t i = 0 ; i < finaljet_list.size() ; ++i) {
        if (i == 7) break;
        _histJetMult->fill(i+1);  // inclusive multiplicity
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const WFinder& wfinder_mu = apply<WFinder>(event, "WFinder_mu");
      if (wfinder_mu.bosons().size() != 1) vetoEvent;

      //const FourMomentum& lepton0 = wfinder_mu.constituentLeptons()[0].momentum();
      //const FourMomentum& neutrino = wfinder_mu.constituentNeutrinos()[0].momentum();
      //double WmT = sqrt( 2 * lepton0.pT() * neutrino.pT() * (1 - cos(deltaPhi(lepton0, neutrino))) );

      const FourMomentum& lepton0 = wfinder_mu.constituentLepton().momentum();
      double WmT = wfinder_mu.mT();

      if (WmT < 50.0*GeV) vetoEvent;
      if (lepton0.abseta() > 2.1 || lepton0.pT() < 25.0*GeV) vetoEvent;

      // Select final jets, ordered by decreasing pT
      vector<FourMomentum> finaljet_list;
      double HT = 0.0;
      const Jets jListAll = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      for (const Jet& j : jListAll) {
        if (j.abseta() < 2.4 && j.pT() > 30.0*GeV && deltaR(lepton0, j) > 0.5) {
          finaljet_list.push_back(j.momentum());
          HT += j.pT();
        }
      }

      // Another jet list, sorted by increasing rapidity
      vector<FourMomentum> jListRap = finaljet_list;
      std::sort(jListRap.begin(), jListRap.end(), cmpMomByRap);

      // Multiplicity exc plot.
      if (finaljet_list.size()<=7) {
        _hist_Mult_exc->fill(finaljet_list.size());
      } else if (finaljet_list.size()>7){
        _hist_Mult_exc->fill(7.);
      }
      // Multiplicity inc plot.
      _fill(_hist_inc_WJetMult, finaljet_list);

      if (finaljet_list.size()>=1) {
        _hist_addJetPt1j->fill(finaljet_list[0].pT());
        _hist_Jeteta1j->fill(fabs(finaljet_list[0].eta()));
        _hist_addHt_1j->fill(HT);
        _hist_dphij1mu_1j->fill( deltaPhi(finaljet_list[0].phi(), lepton0.phi()) );
        _hist_MeanNJht_1j->fill( HT, finaljet_list.size());
      }

      if (finaljet_list.size()>=2) {
        _hist_addJetPt2j->fill(finaljet_list[1].pT());
        _hist_Jeteta2j->fill(fabs(finaljet_list[1].eta()));
        _hist_addHt_2j->fill(HT);

        _hist_dyj1j2_2j   ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()));
        _hist_dyjFjB_2j   ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()));
        _hist_dphij1j2_2j ->fill( deltaPhi(finaljet_list[0].phi(), finaljet_list[1].phi()));
        _hist_dphijFjB_2j ->fill( deltaPhi(jListRap[0].phi(), jListRap[jListRap.size()-1].phi()) );

        _hist_dijetM_2j   ->fill( (add(finaljet_list[0], finaljet_list[1])).mass());
        _hist_diJetPt_2j  ->fill( (add(finaljet_list[0], finaljet_list[1])).pT());
        _hist_dRj1j2_2j   ->fill( deltaR(finaljet_list[0].rapidity(), finaljet_list[0].phi(), finaljet_list[1].rapidity(), finaljet_list[1].phi()));

        _hist_dphij2mu_2j ->fill( deltaPhi(finaljet_list[1].phi(), lepton0.phi()) );

        _hist_MeanNJht_2j->fill( HT, finaljet_list.size());
        _hist_MeanNJdyj1j2_2j->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()), finaljet_list.size());
        _hist_MeanNJdyjFjB_2j->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()), finaljet_list.size());
      }

      if (finaljet_list.size()>=3) {
        _hist_addJetPt3j->fill(finaljet_list[2].pT());
        _hist_Jeteta3j->fill(fabs(finaljet_list[2].eta()));
        _hist_addHt_3j->fill(HT);

        _hist_dyj1j2_3j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()));
        _hist_dyj1j3_3j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[2].rapidity()));
        _hist_dyj2j3_3j     ->fill( fabs(finaljet_list[1].rapidity() - finaljet_list[2].rapidity()));
        _hist_dyjFjB_3j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()));

        _hist_dijetM_3j  ->fill( (add(finaljet_list[0], finaljet_list[1])).mass());
        _hist_diJetPt_3j ->fill( (add(finaljet_list[0], finaljet_list[1])).pT());

        _hist_dphij3mu_3j->fill( deltaPhi(finaljet_list[2].phi(), lepton0.phi()) );
      }

      if (finaljet_list.size()>=4) {
        _hist_addJetPt4j->fill(finaljet_list[3].pT());
        _hist_Jeteta4j->fill(fabs(finaljet_list[3].eta()));
        _hist_addHt_4j->fill(HT);

        _hist_dyj1j2_4j     ->fill( fabs(finaljet_list[0].rapidity() - finaljet_list[1].rapidity()));
        _hist_dyjFjB_4j     ->fill( fabs(jListRap[0].rapidity() - jListRap[jListRap.size()-1].rapidity()));

        _hist_dijetM_4j  ->fill( (add(finaljet_list[0], finaljet_list[1])).mass());
        _hist_diJetPt_4j ->fill( (add(finaljet_list[0], finaljet_list[1])).pT());
        _hist_dphij4mu_4j->fill( deltaPhi(finaljet_list[3].phi(), lepton0.phi()) );
      }

    } //void loop


    /// Normalise histograms etc., after the run
    void finalize() {

      const double crossec = !std::isnan(crossSectionPerEvent()) ? crossSection() : 36703*picobarn;
      if (std::isnan(crossSectionPerEvent())){
        MSG_INFO("No valid cross-section given, using NNLO xsec calculated by FEWZ " << crossec/picobarn << " pb");
      }

      scale(_hist_Mult_exc, crossec/picobarn/sumOfWeights());
      scale(_hist_inc_WJetMult, crossec/picobarn/sumOfWeights());

      scale(_hist_addJetPt1j, crossec/picobarn/sumOfWeights());
      scale(_hist_addJetPt2j, crossec/picobarn/sumOfWeights());
      scale(_hist_addJetPt3j, crossec/picobarn/sumOfWeights());
      scale(_hist_addJetPt4j, crossec/picobarn/sumOfWeights());

      scale(_hist_Jeteta1j, crossec/picobarn/sumOfWeights());
      scale(_hist_Jeteta2j, crossec/picobarn/sumOfWeights());
      scale(_hist_Jeteta3j, crossec/picobarn/sumOfWeights());
      scale(_hist_Jeteta4j, crossec/picobarn/sumOfWeights());

      scale(_hist_addHt_1j, crossec/picobarn/sumOfWeights());
      scale(_hist_addHt_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_addHt_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_addHt_4j, crossec/picobarn/sumOfWeights());

      //-------------------------------------
      scale(_hist_dyj1j2_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dyj1j2_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dyj1j2_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dyjFjB_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dyjFjB_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dyjFjB_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dyj1j3_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dyj2j3_3j, crossec/picobarn/sumOfWeights());

      scale(_hist_dphij1j2_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphijFjB_2j, crossec/picobarn/sumOfWeights());

      scale(_hist_dRj1j2_2j, crossec/picobarn/sumOfWeights());

      scale(_hist_dijetM_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dijetM_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dijetM_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_diJetPt_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_diJetPt_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_diJetPt_4j, crossec/picobarn/sumOfWeights());

      scale(_hist_dphij1mu_1j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij2mu_2j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij3mu_3j, crossec/picobarn/sumOfWeights());
      scale(_hist_dphij4mu_4j, crossec/picobarn/sumOfWeights());

    }

    //@}

  private:

    /// @name Histograms
    //@{

    Histo1DPtr _hist_inc_WJetMult;
    Histo1DPtr _hist_Mult_exc;

    Histo1DPtr _hist_addJetPt1j;
    Histo1DPtr _hist_addJetPt2j;
    Histo1DPtr _hist_addJetPt3j;
    Histo1DPtr _hist_addJetPt4j;

    Histo1DPtr _hist_Jeteta1j;
    Histo1DPtr _hist_Jeteta2j;
    Histo1DPtr _hist_Jeteta3j;
    Histo1DPtr _hist_Jeteta4j;

    Histo1DPtr _hist_addHt_1j;
    Histo1DPtr _hist_addHt_2j;
    Histo1DPtr _hist_addHt_3j;
    Histo1DPtr _hist_addHt_4j;

    //-------------------------------------
    Histo1DPtr _hist_dyj1j2_2j;
    Histo1DPtr _hist_dyj1j2_3j;
    Histo1DPtr _hist_dyj1j2_4j;

    Histo1DPtr _hist_dyjFjB_2j;
    Histo1DPtr _hist_dyjFjB_3j;
    Histo1DPtr _hist_dyjFjB_4j;

    Histo1DPtr _hist_dyj1j3_3j;
    Histo1DPtr _hist_dyj2j3_3j;

    Histo1DPtr _hist_dphij1j2_2j;
    Histo1DPtr _hist_dphijFjB_2j;

    Histo1DPtr _hist_dRj1j2_2j;

    Histo1DPtr _hist_dijetM_2j;
    Histo1DPtr _hist_dijetM_3j;
    Histo1DPtr _hist_dijetM_4j;

    Histo1DPtr _hist_diJetPt_2j;
    Histo1DPtr _hist_diJetPt_3j;
    Histo1DPtr _hist_diJetPt_4j;

    Histo1DPtr _hist_dphij1mu_1j;
    Histo1DPtr _hist_dphij2mu_2j;
    Histo1DPtr _hist_dphij3mu_3j;
    Histo1DPtr _hist_dphij4mu_4j;

    //-------------------------------------
    Profile1DPtr _hist_MeanNJht_1j;
    Profile1DPtr _hist_MeanNJht_2j;
    Profile1DPtr _hist_MeanNJdyj1j2_2j;
    Profile1DPtr _hist_MeanNJdyjFjB_2j;

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1491953);


}
