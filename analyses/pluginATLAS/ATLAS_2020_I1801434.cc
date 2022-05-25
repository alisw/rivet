#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {


  /// @brief All-hadronic ttbar cross-sections at 13 TeV
  class ATLAS_2020_I1801434 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2020_I1801434);


    void init() {

      Cut eta_full =  Cuts::abseta < 5.0;
      Cut lep_cuts = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;

      FinalState fs(eta_full);
      FinalState fs_neutrino;

      const FinalState all_photons(eta_full && Cuts::abspid == PID::PHOTON);
      PromptFinalState photons(all_photons);
      photons.acceptTauDecays(false);
      declare(photons, "photons");


      PromptFinalState electrons(eta_full && Cuts::abspid == PID::ELECTRON);
      electrons.acceptTauDecays(true);
      declare(electrons, "electrons");

      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true, true);
      declare(dressedelectrons, "dressedelectrons");

      DressedLeptons ewdressedelectrons(all_photons, electrons, 0.1, eta_full, true, true);
      declare(ewdressedelectrons, "ewdressedelectrons");

      PromptFinalState muons(eta_full && Cuts::abspid == PID::MUON);
      muons.acceptTauDecays(true);
      declare(muons, "muons");

      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true, true);
      declare(dressedmuons, "dressedmuons");

      DressedLeptons ewdressedmuons(all_photons, muons, 0.1, eta_full, true, true);
      declare(ewdressedmuons, "ewdressedmuons");

      PromptFinalState taus(eta_full && Cuts::abspid == PID::TAU);
      declare(taus, "taus");

      VetoedFinalState vfs(fs);
      InvisibleFinalState prompt_invis(true, true); // require promptness & allow from prompt tau decays
      vfs.addVetoOnThisFinalState(dressedelectrons);
      vfs.addVetoOnThisFinalState(dressedmuons);
      vfs.addVetoOnThisFinalState(prompt_invis);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      declare(jets, "jets");


      /*1*/std::vector<double> jets_n_2D_bins={5.5,6.5,7.5,8.5,9.5};
      /*2*/std::vector<double> mtt_2D_bins={0.0,620.0,835.0,1050.0,3000.0};
      /*3*/std::vector<double> pttop2_2D_bins={0.0, 175.0, 275.0,385.0,1000.0};
      /*4*/std::vector<double> mtt0_2D_bins={0.0,645.0,795.0,1080.0,3000.0};

      book_hist("DR_e1j1",4);
      book_hist("abs_t1_y_1",8);
      book_hist("tt_m",12);
      book_hist("abs_t2_y_1",16);
      book_hist("abs_tt_y",20);
      book_hist("t1_pt",24);
      book_hist("t2_pt",28);
      book_hist("tt_pt",32);
      book_hist("jets_n", 36);
      book_hist("DeltaPhi_1",40);
      book_hist("absPout",44);
      book_hist("absPcross_1",48);
      book_hist("Ztt",52);
      book_hist("HTtt",56);
      book_hist("abs_y_boost",60);
      book_hist("Chitt",64);
      book_hist("RWt1_1",68);
      book_hist("RWt2",72);
      book_hist("RWb1",76);
      book_hist("RWb2",80);
      book_hist("DR_e1tc",84);
      book_hist("DR_e2tc",88);
      book_hist("DR_e3tc",92);
      book_hist("Rpt_e1t1",96);
      book_hist("Rpt_e2t1",100);
      book_hist("Rpt_e3t1",104);
      book_hist("Rpt_tte1",108);
      book_hist("Rpt_e1j1",112);
      book_hist("Rpt_e2j1",116);
      book_hist("Rpt_e3j1",120);
      book_hist("DR_e2e1",124);
      book_hist("DR_e3e1",128);
      book_hist("Rpt_e2e1",132);
      book_hist("Rpt_e3e1",136);


      //--2D--////////////////////////


      book2D("t1_pt_jet_n_multi", jets_n_2D_bins ,153);
      book2D("t1_pt_jet_n_multi_norm", jets_n_2D_bins ,139);
      book2D("t2_pt_jet_n_multi", jets_n_2D_bins , 181);
      book2D("t2_pt_jet_n_multi_norm", jets_n_2D_bins ,167);
      book2D("tt_pt_jet_n_multi", jets_n_2D_bins ,209);
      book2D("tt_pt_jet_n_multi_norm", jets_n_2D_bins ,195);
      book2D("absPout_jet_n_multi", jets_n_2D_bins ,237);
      book2D("absPout_jet_n_multi_norm", jets_n_2D_bins ,223);
      book2D("DeltaPhi_jet_n_multi", jets_n_2D_bins , 265);
      book2D("DeltaPhi_jet_n_multi_norm", jets_n_2D_bins , 251);
      book2D("absPcross_jet_n_multi", jets_n_2D_bins ,293);
      book2D("absPcross_jet_n_multi_norm", jets_n_2D_bins , 279);
      book2D("t2_pt_m_multi", mtt_2D_bins ,321);
      book2D("t2_pt_m_multi_norm", mtt_2D_bins ,307);
      book2D("tt_pt_m_multi", mtt_2D_bins ,349);
      book2D("tt_pt_m_multi_norm", mtt_2D_bins ,335);
      book2D("abs_tt_y_m_multi", mtt_2D_bins ,377);
      book2D("abs_tt_y_m_multi_norm", mtt_2D_bins ,363);
      book2D("t1_pt_t2_pt_multi", pttop2_2D_bins ,405);
      book2D("t1_pt_t2_pt_multi_norm", pttop2_2D_bins ,391);
      book2D("t1_pt_m_multi_y0", mtt0_2D_bins ,433);
      book2D("t1_pt_m_multi_y0_norm", mtt0_2D_bins ,419);

    }


    void analyze(const Event& event) {

      vector<DressedLepton> elecs = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      Particles taus = apply<PromptFinalState>(event, "taus").particlesByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(elecs, jets, 0.4);

      Jets bjets, lightjets;
      for (const Jet& jet: jets){
        bool isBjet = jet.bTagged(Cuts::pT > 5*GeV);
        if (isBjet)  bjets +=jet;
        else         lightjets += jet;
      }

      // Start of the Selection

      //No leptons
      if (elecs.size()) vetoEvent; //No electrons
      if (muons.size()) vetoEvent; //No muons
      if (taus.size())  vetoEvent; //No taus

      //At least 6 jets with pt > 55 GeV
      if (filter_select(jets, Cuts::pT > 55*GeV).size() < 6)  vetoEvent;

      //Exactly 2 bjets
      if ( bjets.size() != 2) vetoEvent;
      //Chi2 calculation
      double minChi2 = 1000.0*TeV;
      const double mWPDG = 80.4*GeV;
      const double sigmaTopSquare = (10.7*GeV)*(10.7*GeV);
      const double sigmaWSquare = (5.9*GeV)*(5.9*GeV);
      int W1j1index = -1, W1j2index = -1;
      int W2j1index = -1, W2j2index = -1;
      int bJet1index = -1, bJet2index = -1;

      for (unsigned int b = 0; b < bjets.size(); ++b) {
        for (unsigned int i = 0; i < (lightjets.size()-1); ++i) {
          for (unsigned int j = i+1; j < (lightjets.size()); ++j) {
            for (unsigned int k = 0; k < (lightjets.size()-1); ++k) {
              for (unsigned int w = k+1; w < lightjets.size(); ++w) {

                FourMomentum W1 = lightjets[i].momentum() + lightjets[j].momentum();
                FourMomentum W2 = lightjets[k].momentum() + lightjets[w].momentum();
                if (lightjets[i].momentum() == lightjets[k].momentum()) continue;
                if (lightjets[i].momentum() == lightjets[w].momentum()) continue;
                if (lightjets[j].momentum() == lightjets[k].momentum()) continue;
                if (lightjets[j].momentum() == lightjets[w].momentum()) continue;

                double wMass1 = W1.mass();
                double wMass2 = W2.mass();

                FourMomentum t1 = bjets[b] + W1;
                FourMomentum t2 = bjets[(b+1)%2] + W2;
                double chi2 = (t1.mass()-t2.mass())*(t1.mass()-t2.mass())/(2*sigmaTopSquare);
                chi2 += (wMass1 - mWPDG)*(wMass1 - mWPDG)/sigmaWSquare;
                chi2 += (wMass2 - mWPDG)*(wMass2 - mWPDG)/sigmaWSquare;


                if (chi2 < minChi2) {
                  minChi2 = chi2;
                  if (t1.pt() > t2.pt()){
                    W1j1index = i;
                    W1j2index = j;
                    W2j1index = k;
                    W2j2index = w;
                    bJet1index = b;
                    bJet2index = (b+1)%2;
                  }
                  else{
                    W1j1index = k;
                    W1j2index = w;
                    W2j1index = i;
                    W2j2index = j;
                    bJet1index = (b+1)%2;
                    bJet2index = b;
                  }
                }
              }
            }
          }
        }
      }


      FourMomentum pw1jet1 = lightjets[W1j1index].momentum();
      FourMomentum pw1jet2 = lightjets[W1j2index].momentum();
      FourMomentum pw2jet1 = lightjets[W2j1index].momentum();
      FourMomentum pw2jet2 = lightjets[W2j2index].momentum();
      FourMomentum bjet1 = bjets[bJet1index].momentum();
      FourMomentum bjet2 = bjets[bJet2index].momentum();

      FourMomentum W1 = pw1jet1 + pw1jet2;
      FourMomentum W2 = pw2jet1 + pw2jet2;
      double drbw1 = deltaR(bjet1,W1);
      double drbw2 = deltaR(bjet2,W2);

      FourMomentum t1 = bjet1 + W1;
      FourMomentum t2 = bjet2 + W2;
      FourMomentum pttbar = t1 + t2;

      //Vector 3
      Vector3 z_versor(0,0,1);
      Vector3 vt1 = t1.vector3();
      Vector3 vt2 = t2.vector3();

      // Variables
      const double HT_ttbar = t1.pt() + t2.pt();
      const double absPout = fabs(vt2.dot((vt1.cross(z_versor))/(vt1.cross(z_versor).mod())));
      size_t jet_multiplicity = jets.size();
      size_t jet_multiplicity_2D = TransformJetMultiplicity(jet_multiplicity);
      const double abs_y1 = t1.absrap();
      const double abs_y2 = t2.absrap();
      const double ystar = (t1.pt() > t2.pt()) ? 0.5 * (t1.rap()-t2.rap()) : 0.5*(t2.rap()-t1.rap());
      const double Chi = exp(2 * fabs(ystar));
      const double Ztt = t2.pt()/t1.pt();
      const double DPhi= deltaPhi(t1,t2);
      const double abs_yboost = fabs(0.5*(t1.rap() +t2.rap()));
      const double RWb1 = W1.pt()/bjet1.pt();
      const double RWb2 = W2.pt()/bjet2.pt();
      const double RWt1 = W1.pt()/t1.pt();
      const double RWt2 = W2.pt()/t2.pt();


      //define extrajets
      vector<int> index_extrajet;
      for (int j = 0; j < int(lightjets.size()); ++j) {
        if (W1j1index != j && W1j2index != j && W2j1index != j && W2j2index != j) index_extrajet.push_back(j);
      }
      double DR_e1j1 = 10000; double DR_e1t1 = 10000; double DR_e1t2 = 10000; double DR_e1tc = 10000;
      double Rpt_e1j1 = -100; double Rpt_e1t1 = -100; double Rpt_tte1 = 100;
      double DR_e2t1 = 10000; double DR_e2t2 = 10000; double DR_e2tc = 10000;
      double Rpt_e2j1 = -100; double Rpt_e2t1 = -100;
      double DR_e3t1 = 10000; double DR_e3t2 = 10000; double DR_e3tc = 10000;
      double Rpt_e3j1 = -100; double Rpt_e3t1 = -100;
      double DR_e2e1 = 10000; double DR_e3e1 = 10000;
      double Rpt_e2e1 = 100; double Rpt_e3e1 = 100;

      if (index_extrajet.size()) {
        DR_e1j1 = deltaR(lightjets[index_extrajet.at(0)], jets[0]);
        DR_e1t1 = deltaR(lightjets[index_extrajet.at(0)], t1);
        DR_e1t2 = deltaR(lightjets[index_extrajet.at(0)], t2);
        Rpt_e1j1 = lightjets[index_extrajet.at(0)].pt()/jets[0].pt();
        Rpt_e1t1 = lightjets[index_extrajet.at(0)].pt()/t1.pt();
        Rpt_tte1 = deltaR(pttbar, lightjets[index_extrajet.at(0)]);
        if (DR_e1t1>DR_e1t2)  DR_e1tc = DR_e1t2;
        else                  DR_e1tc = DR_e1t1;
      }

      if (index_extrajet.size() > 1) {

        DR_e2t1 = deltaR(lightjets[index_extrajet.at(1)], t1);
        DR_e2t2 = deltaR(lightjets[index_extrajet.at(1)], t2);
        Rpt_e2j1 = lightjets[index_extrajet.at(1)].pt()/jets[0].pt();
        Rpt_e2t1 = lightjets[index_extrajet.at(1)].pt()/t1.pt();
        if(DR_e2t1 > DR_e2t2)  DR_e2tc = DR_e2t2;
        else                   DR_e2tc = DR_e2t1;
        Rpt_e2e1 = lightjets[index_extrajet.at(1)].pt()/lightjets[index_extrajet.at(0)].pt();
        DR_e2e1 = deltaR(lightjets[index_extrajet.at(1)],lightjets[index_extrajet.at(0)]);

      }


      if (index_extrajet.size() > 2) {

        DR_e3t1 = deltaR(lightjets[index_extrajet.at(2)], t1);
        DR_e3t2 = deltaR(lightjets[index_extrajet.at(2)], t2);
        Rpt_e3j1 = lightjets[index_extrajet.at(2)].pt()/jets[0].pt();
        Rpt_e3t1 = lightjets[index_extrajet.at(2)].pt()/t1.pt();
        if (DR_e3t1>DR_e3t2) DR_e3tc  = DR_e3t2;
        else                 DR_e3tc  = DR_e3t1;
        Rpt_e3e1 = lightjets[index_extrajet.at(2)].pt()/lightjets[index_extrajet.at(0)].pt();
        DR_e3e1 = deltaR(lightjets[index_extrajet.at(2)],lightjets[index_extrajet.at(0)]);
      }

      //Cut on minChi2
      double absPcross=fabs(p_cross(lightjets[W1j1index], lightjets[W1j2index],
                                    lightjets[W2j1index], lightjets[W2j2index],
                                    bjets[bJet1index],    bjets[bJet2index]));

      if (minChi2 > 10) vetoEvent;
      //Cut on dRbb
      if (deltaR(bjet1,bjet2) < 2.0) vetoEvent;

      //Cut on max dR(b,W)
      if (max(drbw1,drbw2) > 2.2) vetoEvent;

      // Cut on masses
      if (t1.mass() < 130*GeV || t1.mass() >= 200*GeV) vetoEvent;
      if (t2.mass() < 130*GeV || t2.mass() >= 200*GeV) vetoEvent;


      _h["t1_pt"]->fill(t1.pt()/GeV);
      _h["t2_pt"]->fill(t2.pt()/GeV);
      _h["tt_pt"]->fill(pttbar.pt()/GeV);
      _h["absPout"]->fill(absPout);
      _h["jets_n"]->fill(jet_multiplicity);
      _h["abs_t1_y_1"]->fill(abs_y1);
      _h["abs_t2_y_1"]->fill(abs_y2);
      _h["abs_tt_y"]->fill(pttbar.absrap());
      _h["tt_m"]->fill(pttbar.mass()/GeV);
      _h["HTtt"]->fill(HT_ttbar/GeV);
      _h["Chitt"]->fill(Chi);
      _h["Ztt"]->fill(Ztt);
      _h["DeltaPhi_1"]->fill(DPhi);
      _h["abs_y_boost"]->fill(abs_yboost);
      _h["absPcross_1"]->fill(absPcross);
      _h["RWb1"]->fill(RWb1);
      _h["RWb2"]->fill(RWb2);
      _h["RWt1_1"]->fill(RWt1);
      _h["RWt2"]->fill(RWt2);
      _h["Rpt_tte1"]->fill(Rpt_tte1);
      _h["Rpt_e1t1"]->fill(Rpt_e1t1);
      _h["DR_e1tc"]->fill(DR_e1tc);
      _h["Rpt_e2t1"]->fill(Rpt_e2t1);
      _h["DR_e2tc"]->fill(DR_e2tc);
      _h["Rpt_e3t1"]->fill(Rpt_e3t1);
      _h["DR_e3tc"]->fill(DR_e3tc);
      _h["Rpt_e1j1"]->fill(Rpt_e1j1);
      _h["Rpt_e2j1"]->fill(Rpt_e2j1);
      _h["Rpt_e3j1"]->fill(Rpt_e3j1);
      _h["Rpt_e2e1"]->fill(Rpt_e2e1);
      _h["Rpt_e3e1"]->fill(Rpt_e3e1);
      _h["DR_e1j1"]->fill(DR_e1j1);
      _h["DR_e2e1"]->fill(DR_e2e1);
      _h["DR_e3e1"]->fill(DR_e3e1);

      _h["t1_pt_norm"]->fill(t1.pt()/GeV);
      _h["t2_pt_norm"]->fill(t2.pt()/GeV);
      _h["tt_pt_norm"]->fill(pttbar.pt()/GeV);
      _h["absPout_norm"]->fill(absPout);
      _h["jets_n_norm"]->fill(jet_multiplicity);
      _h["abs_t1_y_1_norm"]->fill(abs_y1);
      _h["abs_t2_y_1_norm"]->fill(abs_y2);
      _h["abs_tt_y_norm"]->fill(pttbar.absrap());
      _h["tt_m_norm"]->fill(pttbar.mass()/GeV);
      _h["HTtt_norm"]->fill(HT_ttbar/GeV);
      _h["Chitt_norm"]->fill(Chi);
      _h["Ztt_norm"]->fill(Ztt);
      _h["DeltaPhi_1_norm"]->fill(DPhi);
      _h["abs_y_boost_norm"]->fill(abs_yboost);
      _h["absPcross_1_norm"]->fill(absPcross);
      _h["RWb1_norm"]->fill(RWb1);
      _h["RWb2_norm"]->fill(RWb2);
      _h["RWt1_1_norm"]->fill(RWt1);
      _h["RWt2_norm"]->fill(RWt2);
      _h["Rpt_tte1_norm"]->fill(Rpt_tte1);
      _h["Rpt_e1t1_norm"]->fill(Rpt_e1t1);
      _h["DR_e1tc_norm"]->fill(DR_e1tc);
      _h["Rpt_e2t1_norm"]->fill(Rpt_e2t1);
      _h["DR_e2tc_norm"]->fill(DR_e2tc);
      _h["Rpt_e3t1_norm"]->fill(Rpt_e3t1);
      _h["DR_e3tc_norm"]->fill(DR_e3tc);
      _h["Rpt_e1j1_norm"]->fill(Rpt_e1j1);
      _h["Rpt_e2j1_norm"]->fill(Rpt_e2j1);
      _h["Rpt_e3j1_norm"]->fill(Rpt_e3j1);
      _h["Rpt_e2e1_norm"]->fill(Rpt_e2e1);
      _h["Rpt_e3e1_norm"]->fill(Rpt_e3e1);
      _h["DR_e1j1_norm"]->fill(DR_e1j1);
      _h["DR_e2e1_norm"]->fill(DR_e2e1);
      _h["DR_e3e1_norm"]->fill(DR_e3e1);


      _h_multi["t1_pt_jet_n_multi"].fill(jet_multiplicity_2D, t1.pt()/GeV);
      _h_multi["t2_pt_jet_n_multi"].fill(jet_multiplicity_2D, t2.pt()/GeV);
      _h_multi["tt_pt_jet_n_multi"].fill(jet_multiplicity_2D, pttbar.pt()/GeV);
      _h_multi["absPout_jet_n_multi"].fill(jet_multiplicity_2D, absPout);
      _h_multi["DeltaPhi_jet_n_multi"].fill(jet_multiplicity_2D, DPhi);
      _h_multi["absPcross_jet_n_multi"].fill(jet_multiplicity_2D, absPcross);
      _h_multi["t2_pt_m_multi"].fill(pttbar.mass()/GeV, t2.pt()/GeV);
      _h_multi["tt_pt_m_multi"].fill(pttbar.mass()/GeV, pttbar.pt()/GeV);
      _h_multi["abs_tt_y_m_multi"].fill(pttbar.mass()/GeV, pttbar.absrap());
      _h_multi["t1_pt_t2_pt_multi"].fill(t2.pt()/GeV, t1.pt()/GeV);
      _h_multi["t1_pt_m_multi_y0"].fill(pttbar.mass()/GeV, t1.pt()/GeV);
      _h_multi["t1_pt_jet_n_multi_norm"].fill(jet_multiplicity_2D, t1.pt()/GeV);
      _h_multi["t2_pt_jet_n_multi_norm"].fill(jet_multiplicity_2D, t2.pt()/GeV);
      _h_multi["tt_pt_jet_n_multi_norm"].fill(jet_multiplicity_2D, pttbar.pt()/GeV);
      _h_multi["absPout_jet_n_multi_norm"].fill(jet_multiplicity_2D, absPout);
      _h_multi["DeltaPhi_jet_n_multi_norm"].fill(jet_multiplicity_2D, DPhi);
      _h_multi["absPcross_jet_n_multi_norm"].fill(jet_multiplicity_2D, absPcross);
      _h_multi["t2_pt_m_multi_norm"].fill(pttbar.mass()/GeV, t2.pt()/GeV);
      _h_multi["tt_pt_m_multi_norm"].fill(pttbar.mass()/GeV, pttbar.pt()/GeV);
      _h_multi["abs_tt_y_m_multi_norm"].fill(pttbar.mass()/GeV, pttbar.absrap());
      _h_multi["t1_pt_t2_pt_multi_norm"].fill(t2.pt()/GeV, t1.pt()/GeV);
      _h_multi["t1_pt_m_multi_y0_norm"].fill(pttbar.mass()/GeV, t1.pt()/GeV);
    }

    void finalize() {

      // Normalize to cross-section
      const double sf = crossSection() / sumOfWeights();
      for (auto& hit : _h) {
        scale(hit.second, sf);
        if (hit.first.find("_norm") != string::npos)  normalize(hit.second, 1.0, false);
      }
      for (auto& hit : _h_multi) {
        if (hit.first.find("_norm") != string::npos) {
          for (Histo1DPtr& hist : hit.second.histos()) { scale(hist, sf); }
          const double norm2D = integral2D(hit.second);
          hit.second.scale(safediv(1.0, norm2D), this);
        }
        else {
          hit.second.scale(sf, this);
        }
      }
    }

  private:

    double integral2D(BinnedHistogram& h_multi) {
      double total_integral = 0;
      for  (Histo1DPtr& h : h_multi.histos()) {
        total_integral += h->integral(false);
      }
      return total_integral;
    }

    void book2D(string name, std::vector<double>& doubleDiff_bins, size_t table){
      for (size_t i = 0; i < doubleDiff_bins.size() - 1; ++i) {
        string nbin = std::to_string(i);
        { Histo1DPtr tmp; _h_multi[name].add(doubleDiff_bins[i], doubleDiff_bins[i+1], book(tmp, table+i, 1, 1)); }
      }
    }

    void book_hist(string name, size_t table) {
      book(_h[name], table, 1, 1);
      book(_h[name+"_norm"], table-2, 1, 1);
    }


    int TransformJetMultiplicity(int jet_n){
      int new_jet_n = -1;
      if (jet_n >= 9)  new_jet_n = 9;
      else             new_jet_n = jet_n;
      return new_jet_n;
    }

    double p_cross(FourMomentum j1, FourMomentum j2, FourMomentum j3, FourMomentum j4, FourMomentum b1, FourMomentum b2) {
      Vector3 vj1 = j1.vector3().unit();
      Vector3 vj2 = j2.vector3().unit();
      Vector3 vj3 = j3.vector3().unit();
      Vector3 vj4 = j4.vector3().unit();
      Vector3 vb1 = b1.vector3().unit();
      Vector3 vb2 = b2.vector3().unit();
      vj1.mod();
      vj2.mod();
      vj3.mod();
      vj4.mod();
      vb1.mod();
      vb2.mod();
      Vector3 vj1j2 = vj1.cross( vj2 );
      Vector3 vj3j4 = vj3.cross( vj4 );
      Vector3 vb1j = vb1.cross( vj1j2 );
      Vector3 vb2j = vb2.cross( vj3j4 );
      Vector3 vcross = vb1j.cross( vb2j );

      return vcross.mod();

    }

    /// @name Objects that are used by the event selection decisions
    map<string, Histo1DPtr> _h;
    map<string, BinnedHistogram> _h_multi;
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2020_I1801434);

}
