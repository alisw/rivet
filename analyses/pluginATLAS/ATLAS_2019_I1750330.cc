// -*- C++ -*
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"


namespace Rivet {

  /// @brief Semileptonic ttbar at 13 TeV
  class ATLAS_2019_I1750330 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1750330);

    void init() {

      _doBoosted = true, _doResolved = true;
      if ( getOption("TYPE") == "BOOSTED" )  _doResolved = false;
      else if ( getOption("TYPE") == "RESOLVED" )  _doBoosted = false;

      Cut eta_full =  (Cuts::abseta < 5.0);
      Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 27*GeV);
      const FinalState fs(eta_full);

      FinalState all_photons(fs, Cuts::abspid == PID::PHOTON);

      PromptFinalState photons(all_photons, false);
      declare(photons, "photons");

      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON, true);
      declare(electrons, "electrons");

      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true);
      declare(dressedelectrons, "dressedelectrons");

      DressedLeptons ewdressedelectrons(all_photons, electrons, 0.1, eta_full, true);
      declare(ewdressedelectrons, "ewdressedelectrons");

      PromptFinalState muons(Cuts::abspid == PID::MUON, true);
      declare(muons, "muons");

      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true);
      declare(dressedmuons, "dressedmuons");

      DressedLeptons ewdressedmuons(all_photons, muons, 0.1, eta_full, true);
      declare(ewdressedmuons, "ewdressedmuons");

      const InvisibleFinalState neutrinos(true, true);

      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(dressedelectrons);
      vfs.addVetoOnThisFinalState(dressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(jets, "boosted_jets");

      VetoedFinalState vfs_res(fs);
      vfs_res.addVetoOnThisFinalState(ewdressedelectrons);
      vfs_res.addVetoOnThisFinalState(ewdressedmuons);
      vfs_res.addVetoOnThisFinalState(neutrinos);
      FastJets jets_res(vfs_res, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(jets_res, "resolved_jets");

      declare(MissingMomentum(), "MissingMomentum");

      // Bins for 2D resolved
      std::vector<double> ttbar_m_2D_bins = {200,400,550,700,1000,2000};
      std::vector<double> top_had_pt_2D_bins = {0,60,120,200,300,1000};
      std::vector<double> ttbar_pt_2D_bins = {0,30,80,190,800};
      std::vector<double> top_had_abs_y_2D_bins = {0,0.7,1.4,2.5};
      std::vector<double> ttbar_abs_y_2D_bins = {0.0,0.4,0.8,1.2,2.5};

      std::vector<double> n_jet_bins = {3.5, 4.5, 5.5, 6.5, 7.5};
      std::vector<double> n_jet_bins_for_ttbar_m = {3.5, 4.5, 5.5, 6.5};
      std::vector<double> n_extrajet_bins = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5};

      //Bins for 2D boosted
      std::vector<double> eta_2D_bins = {0,1,2};
      std::vector<double> etattbar_2D_bins = {0,60,120,200,300,1000};
      std::vector<double> pttbar_2D_bins = {0.0, 40.0, 150.0, 1000.0};
      std::vector<double> mtt_2D_bins = {490.0, 1160, 3000.0};
      std::vector<double> eta_external_2D_bins = {0.0, 0.65, 1.3, 2.0};
      std::vector<double> ptt_external_mtt_2D_bins = {0.0, 40.0, 150.0, 1000.0};
      std::vector<double> Htt_external_2D_bins = {350.0, 780.0, 2500.0};
      std::vector<double> eta_external_ptt_2D_bins = {0.0, 0.65, 2.0};

      std::vector<double> n_jet_pttop_bins = {-0.5,1.5,2.5,3.5};
      std::vector<double> n_jet_ptttbar_bins = {-0.5,1.5,3.5};
      std::vector<double> n_jet_Pout_bins = {-0.5,1.5,3.5};
      std::vector<double> n_jet_mtt_bins = {-0.5,0.5,1.5,2.5};

      //Resolved histograms (digits correspond to "Table ID" from HepData)
      book2D("ttbar_m_top_had_pt_multi_norm", ttbar_m_2D_bins,54);
      book2D("ttbar_m_top_had_pt_multi", ttbar_m_2D_bins, 74);

      book2D("ttbar_m_ttbar_pt_multi_norm", ttbar_m_2D_bins, 94);
      book2D("ttbar_m_ttbar_pt_multi", ttbar_m_2D_bins, 114);

      book2D("top_had_pt_absPout_multi_norm", top_had_pt_2D_bins, 134);
      book2D("top_had_pt_absPout_multi", top_had_pt_2D_bins,154);

      book2D("top_had_pt_jet_n_multi_norm", n_jet_bins,174);
      book2D("top_had_pt_jet_n_multi", n_jet_bins, 188 );

      book2D("ttbar_m_jet_n_multi_norm", n_jet_bins_for_ttbar_m,202);
      book2D("ttbar_m_jet_n_multi", n_jet_bins_for_ttbar_m,211);

      book2D("ttbar_pt_jet_n_multi_norm", n_jet_bins, 220 );
      book2D("ttbar_pt_jet_n_multi", n_jet_bins, 234 );

      book2D("absPout_jet_n_multi_norm", n_jet_bins, 248 );
      book2D("absPout_jet_n_multi", n_jet_bins, 262 );

      book2D("deltaPhi_tt_jet_n_multi_norm", n_jet_bins, 276 );
      book2D("deltaPhi_tt_jet_n_multi", n_jet_bins, 290 );

      book2D("HT_tt_jet_n_multi_norm", n_jet_bins, 304 );
      book2D("HT_tt_jet_n_multi", n_jet_bins, 318 );

      book2D("top_had_abs_y_jet_n_multi_norm", n_jet_bins, 332 );
      book2D("top_had_abs_y_jet_n_multi", n_jet_bins, 346 );

      book2D("ttbar_abs_y_jet_n_multi_norm", n_jet_bins, 360 );
      book2D("ttbar_abs_y_jet_n_multi", n_jet_bins, 374 );

      book2D("chi_tt_jet_n_multi_norm", n_jet_bins, 388 );
      book2D("chi_tt_jet_n_multi", n_jet_bins, 402 );

      book2D("top_had_abs_y_top_had_pt_multi_norm", top_had_abs_y_2D_bins,416 );
      book2D("top_had_abs_y_top_had_pt_multi", top_had_abs_y_2D_bins, 425);

      book2D("ttbar_abs_y_ttbar_pt_multi_norm", ttbar_abs_y_2D_bins, 434);
      book2D("ttbar_abs_y_ttbar_pt_multi", ttbar_abs_y_2D_bins, 448);

      book2D("ttbar_abs_y_ttbar_m_multi_norm", ttbar_abs_y_2D_bins, 462);
      book2D("ttbar_abs_y_ttbar_m_multi", ttbar_abs_y_2D_bins, 476);

      book2D("ttbar_pt_top_had_pt_multi_norm", ttbar_pt_2D_bins, 490);
      book2D("ttbar_pt_top_had_pt_multi", ttbar_pt_2D_bins, 504);

      book_hist("top_had_pt",1);
      book_hist("top_had_abs_y_fine",5);
      book_hist("leading_top_pt",9);
      book_hist("subleading_top_pt",13);
      book_hist("ttbar_m",17);
      book_hist("ttbar_pt",21);
      book_hist("absPout",25);
      book_hist("deltaPhi_tt",29);
      book_hist("HT_tt",33);
      book_hist("extrajet_n",37);
      book_hist("ttbar_abs_y_fine",41);
      book_hist("abs_y_boost",45);
      book_hist("chi_tt",49);

      //Boosted histograms (digits correspond to "Table ID" from HepData)
      book2D("boosted_rc_pttop_etatop_multi", eta_2D_bins, 922);
      book2D("boosted_rc_pttop_etattbar_multi", eta_2D_bins, 912);
      book2D("boosted_rc_pttop_ptttbar_multi", pttbar_2D_bins, 898);
      book2D("boosted_rc_pttop_mttbar_multi", mtt_2D_bins, 932);
      book2D("boosted_rc_mttbar_etattbar_multi", eta_external_2D_bins, 974);
      book2D("boosted_rc_mttbar_ptttbar_multi", ptt_external_mtt_2D_bins, 956);
      book2D("boosted_rc_mttbar_HT_multi", Htt_external_2D_bins, 942);
      book2D("boosted_rc_pttop_extrajet_multi", n_jet_pttop_bins, 992);
      book2D("boosted_rc_ptttbar_extrajet_multi", n_jet_ptttbar_bins, 1006);
      book2D("boosted_rc_mttbar_extrajet_multi", n_jet_mtt_bins, 1020);

      book2D("boosted_rc_pttop_etatop_multi_norm", eta_2D_bins, 917);
      book2D("boosted_rc_pttop_etattbar_multi_norm", eta_2D_bins, 907);
      book2D("boosted_rc_pttop_ptttbar_multi_norm", pttbar_2D_bins, 889);
      book2D("boosted_rc_pttop_mttbar_multi_norm", mtt_2D_bins, 927);
      book2D("boosted_rc_mttbar_etattbar_multi_norm", eta_external_2D_bins, 965);
      book2D("boosted_rc_mttbar_ptttbar_multi_norm", ptt_external_mtt_2D_bins, 947);
      book2D("boosted_rc_mttbar_HT_multi_norm", Htt_external_2D_bins, 937);
      book2D("boosted_rc_pttop_extrajet_multi_norm", n_jet_pttop_bins, 983);
      book2D("boosted_rc_ptttbar_extrajet_multi_norm", n_jet_ptttbar_bins, 1001);
      book2D("boosted_rc_mttbar_extrajet_multi_norm", n_jet_mtt_bins, 1011);

      book_hist("hadTop_boosted_rc_pt",840);
      book_hist("hadTop_boosted_rc_y",844);
      book_hist("LeadingTop_boosted_rc_pt",848);
      book_hist("SubLeadingTop_boosted_rc_pt",852);
      book_hist("boosted_rc_Pout_lep",872);
      book_hist("boosted_rc_chi_tt",868);
      book_hist("boosted_rc_HT",876);
      book_hist("hadTop_boosted_rc_subjets",884);
      book_hist("boosted_rc_extrajet",880);
      book_hist("ttbar_boosted_rc_m",864);
      book_hist("ttbar_boosted_rc_pt",856);
      book_hist("ttbar_boosted_rc_Rapidity",860);

    }


    void analyze(const Event& event) {
      if (_doResolved)  Resolved_selection(event);
      if (_doBoosted)   Boosted_selection(event);
    }


    void Resolved_selection(const Event& event) {

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons     = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      const Jets& jets = apply<FastJets>(event, "resolved_jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      FourMomentum met = apply<MissingMomentum>(event, "MissingMomentum").missingMomentum();

      Jets bjets, lightjets;

      // OVERLAP REMOVAL
      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(electrons, jets, 0.4);

      // b-tagging
      // If there are more than 2 b-tagged jets, the extra b-tagged jets will be treat as light jets
      for (const Jet& jet : jets) {
        bool b_tagged = jet.bTagged(Cuts::pT > 5*GeV);
        if ( b_tagged && bjets.size() < 2)  bjets +=jet;
        else lightjets += jet;
      }

      bool single_electron = electrons.size() == 1 && muons.empty();
      bool single_muon     = muons.size() == 1 && electrons.empty();

      DressedLepton *lepton = NULL;
      if (single_electron)   lepton = &electrons[0];
      else if (single_muon)  lepton = &muons[0];

      if (!single_electron && !single_muon) vetoEvent;
      bool num_b_tagged_jets = (bjets.size() == 2);
      if (!num_b_tagged_jets) vetoEvent;

      if (lightjets.size() < 2) vetoEvent;

      FourMomentum pbjet1; //Momentum of bjet1
      FourMomentum pbjet2; //Momentum of bjet
      if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
        pbjet1 = bjets[0].momentum();
        pbjet2 = bjets[1].momentum();
      } else {
        pbjet1 = bjets[1].momentum();
        pbjet2 = bjets[0].momentum();
      }

      double bestWmass = 1000.0*TeV;
      double mWPDG = 80.399*GeV;
      int Wj1index = -1, Wj2index = -1;
      for (unsigned int i = 0; i < (lightjets.size() - 1); ++i) {
        for (unsigned int j = i + 1; j < lightjets.size(); ++j) {
          double wmass = (lightjets[i].momentum() + lightjets[j].momentum()).mass();
          if (fabs(wmass - mWPDG) < fabs(bestWmass - mWPDG)) {
            bestWmass = wmass;
            Wj1index = i;
            Wj2index = j;
          }
        }
      }

      FourMomentum pjet1 = lightjets[Wj1index].momentum();
      FourMomentum pjet2 = lightjets[Wj2index].momentum();

      // compute hadronic W boson
      FourMomentum pWhadron = pjet1 + pjet2;
      double pz = computeneutrinoz(lepton->momentum(), met);
      FourMomentum ppseudoneutrino( sqrt(sqr(met.px()) + sqr(met.py()) + sqr(pz)), met.px(), met.py(), pz);

      //compute leptonic, hadronic, combined pseudo-top
      FourMomentum ppseudotoplepton = lepton->momentum() + ppseudoneutrino + pbjet1;
      FourMomentum ppseudotophadron = pbjet2 + pWhadron;
      FourMomentum pttbar = ppseudotoplepton + ppseudotophadron;

      Vector3 z_versor(0,0,1);
      Vector3 vpseudotophadron = ppseudotophadron.vector3();
      Vector3 vpseudotoplepton = ppseudotoplepton.vector3();

      // Variables
      double ystar = (ppseudotophadron.pt() > ppseudotoplepton.pt()) ? 0.5 * (ppseudotophadron.rap()-ppseudotoplepton.rap()) : 0.5*(ppseudotoplepton.rap()-ppseudotophadron.rap());
      double chi_ttbar = exp(2 * fabs(ystar));
      double deltaPhi_ttbar = deltaPhi(ppseudotoplepton,ppseudotophadron);
      double HT_ttbar = ppseudotophadron.pt() + ppseudotoplepton.pt();
      double Yboost = 0.5 * (ppseudotophadron.rapidity() + ppseudotoplepton.rapidity());
      double Pout = vpseudotophadron.dot((vpseudotoplepton.cross(z_versor))/(vpseudotoplepton.cross(z_versor).mod()));
      double absPout = fabs(Pout);
      double Leading_top_pt = (ppseudotophadron.pt() > ppseudotoplepton.pt()) ? ppseudotophadron.pt() : ppseudotoplepton.pt();
      double Subleading_top_pt = (ppseudotophadron.pt() > ppseudotoplepton.pt()) ? ppseudotoplepton.pt() : ppseudotophadron.pt();
      int jet_multiplicity = jets.size();
      int extrajet_n = jet_multiplicity - 4;
      int new_jet_multi = TransformJetMultiplicity(jet_multiplicity);
      int new_jet_multi_for_ttbar_m = TransformJetMultiplicity_for_ttbar_m(jet_multiplicity);
      int new_extrajet_multi =TransformExtrajetMultiplicity(extrajet_n);

      _h_multi["top_had_pt_Pout_multi"].fill(ppseudotophadron.pt()/GeV, Pout);
      _h_multi["top_had_pt_absPout_multi"].fill(ppseudotophadron.pt()/GeV, absPout);
      _h_multi["ttbar_m_top_had_pt_multi"].fill(pttbar.mass()/GeV, ppseudotophadron.pt()/GeV);
      _h_multi["ttbar_m_ttbar_pt_multi"].fill(pttbar.mass()/GeV, pttbar.pt()/GeV);
      _h_multi["ttbar_pt_top_had_pt_multi"].fill(pttbar.pt()/GeV, ppseudotophadron.pt()/GeV);
      _h_multi["ttbar_abs_y_ttbar_pt_multi"].fill(pttbar.absrap(), pttbar.pt()/GeV);
      _h_multi["ttbar_abs_y_ttbar_m_multi"].fill(pttbar.absrap(), pttbar.mass()/GeV);
      _h_multi["top_had_abs_y_top_had_pt_multi"].fill(ppseudotophadron.absrap(), ppseudotophadron.pt()/GeV);

      _h_multi["ttbar_pt_jet_n_multi"].fill(new_jet_multi, pttbar.pt()/GeV);
      _h_multi["ttbar_m_jet_n_multi"].fill(new_jet_multi_for_ttbar_m, pttbar.mass()/GeV);
      _h_multi["chi_tt_jet_n_multi"].fill(new_jet_multi, chi_ttbar);
      _h_multi["absPout_jet_n_multi"].fill(new_jet_multi, absPout);
      _h_multi["deltaPhi_tt_jet_n_multi"].fill(new_jet_multi, deltaPhi_ttbar);
      _h_multi["HT_tt_jet_n_multi"].fill(new_jet_multi, HT_ttbar/GeV);
      _h_multi["top_had_pt_jet_n_multi"].fill(new_jet_multi, ppseudotophadron.pt()/GeV);
      _h_multi["top_had_abs_y_jet_n_multi"].fill(new_jet_multi, ppseudotophadron.absrap());
      _h_multi["ttbar_abs_y_jet_n_multi"].fill(new_jet_multi, pttbar.absrap());

      _h_multi["top_had_pt_Pout_multi_norm"].fill(ppseudotophadron.pt()/GeV, Pout);
      _h_multi["top_had_pt_absPout_multi_norm"].fill(ppseudotophadron.pt()/GeV, absPout);
      _h_multi["ttbar_m_top_had_pt_multi_norm"].fill(pttbar.mass()/GeV, ppseudotophadron.pt()/GeV);
      _h_multi["ttbar_m_ttbar_pt_multi_norm"].fill(pttbar.mass()/GeV, pttbar.pt()/GeV);
      _h_multi["ttbar_pt_top_had_pt_multi_norm"].fill(pttbar.pt()/GeV, ppseudotophadron.pt()/GeV);
      _h_multi["ttbar_abs_y_ttbar_pt_multi_norm"].fill(pttbar.absrap(), pttbar.pt()/GeV);
      _h_multi["ttbar_abs_y_ttbar_m_multi_norm"].fill(pttbar.absrap(), pttbar.mass()/GeV);
      _h_multi["top_had_abs_y_top_had_pt_multi_norm"].fill(ppseudotophadron.absrap(), ppseudotophadron.pt()/GeV);

      _h_multi["ttbar_pt_jet_n_multi_norm"].fill(new_jet_multi, pttbar.pt());
      _h_multi["ttbar_m_jet_n_multi_norm"].fill(new_jet_multi_for_ttbar_m, pttbar.mass());
      _h_multi["chi_tt_jet_n_multi_norm"].fill(new_jet_multi, chi_ttbar);
      _h_multi["absPout_jet_n_multi_norm"].fill(new_jet_multi, absPout);
      _h_multi["deltaPhi_tt_jet_n_multi_norm"].fill(new_jet_multi, deltaPhi_ttbar);
      _h_multi["HT_tt_jet_n_multi_norm"].fill(new_jet_multi, HT_ttbar/GeV);
      _h_multi["top_had_pt_jet_n_multi_norm"].fill(new_jet_multi, ppseudotophadron.pt()/GeV);
      _h_multi["top_had_abs_y_jet_n_multi_norm"].fill(new_jet_multi, ppseudotophadron.absrap());
      _h_multi["ttbar_abs_y_jet_n_multi_norm"].fill(new_jet_multi, pttbar.absrap());

      _h["chi_tt"]->fill(chi_ttbar);
      _h["deltaPhi_tt"]->fill(deltaPhi_ttbar);
      _h["HT_tt"]->fill(HT_ttbar/GeV);
      _h["absPout"]->fill(absPout);
      _h["abs_y_boost"]->fill(fabs(Yboost));
      _h["top_had_pt"]->fill(ppseudotophadron.pt()/GeV);
      _h["top_had_abs_y_fine"]->fill(ppseudotophadron.absrap());
      _h["ttbar_pt"]->fill(pttbar.pt()/GeV);
      _h["ttbar_m"]->fill(pttbar.mass()/GeV);
      _h["ttbar_abs_y_fine"]->fill(pttbar.absrap());
      _h["leading_top_pt"]->fill(Leading_top_pt/GeV);
      _h["subleading_top_pt"]->fill(Subleading_top_pt/GeV);
      _h["extrajet_n"]->fill(new_extrajet_multi+1);

      _h["chi_tt_norm"]->fill(chi_ttbar);
      _h["deltaPhi_tt_norm"]->fill(deltaPhi_ttbar);
      _h["HT_tt_norm"]->fill(HT_ttbar/GeV);
      _h["absPout_norm"]->fill(absPout);
      _h["abs_y_boost_norm"]->fill(fabs(Yboost));
      _h["top_had_pt_norm"]->fill(ppseudotophadron.pt()/GeV);
      _h["top_had_abs_y_fine_norm"]->fill(ppseudotophadron.absrap());
      _h["ttbar_pt_norm"]->fill(pttbar.pt()/GeV);
      _h["ttbar_m_norm"]->fill(pttbar.mass()/GeV);
      _h["ttbar_abs_y_fine_norm"]->fill(pttbar.absrap());
      _h["leading_top_pt_norm"]->fill(Leading_top_pt/GeV);
      _h["subleading_top_pt_norm"]->fill(Subleading_top_pt/GeV);
      _h["extrajet_n_norm"]->fill(new_extrajet_multi+1);
    }

    void Boosted_selection(const Event& event) {

      //Projections
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      const Jets& jets = apply<FastJets>(event, "boosted_jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta <= 2.5);
      const FourMomentum& met = apply<MissingMomentum>(event, "MissingMomentum").missingMomentum();

      if (jets.size() < 2)  vetoEvent;
      PseudoJets smallRjets;
      for (const Jet& jet : jets) {
        smallRjets += jet.pseudojet();
        bool b_tagged = jet.bTagged(Cuts::pT > 5*GeV);
        smallRjets[smallRjets.size()-1].set_user_index(b_tagged); // cheeky, but works
      }

      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(electrons, jets, 0.4);

      fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0), fastjet::SelectorPtFractionMin(0.05));
      fastjet::ClusterSequence antikt_cs(smallRjets, fastjet::JetDefinition(fastjet::antikt_algorithm, 1.0));
      PseudoJets reclustered_jets = antikt_cs.inclusive_jets();

      // trim the jets
      Jets TrimmedJets;
      for (PseudoJet pjet : reclustered_jets) {
        PseudoJet ptrim = trimmer(pjet);
        if (ptrim.perp() < 350*GeV)  continue;
        if (fabs(ptrim.eta()) > 2.0) continue;
        bool bTagged  = false;
        Particles constituents;
        for (const PseudoJet& c : ptrim.constituents()) {
          // we only care about the number of subjets, so
          // fine to treat as Particles with dummy PID
          constituents += Particle(0, momentum(c));
          bTagged |= c.user_index();
        }
        ptrim.set_user_index(bTagged);
        TrimmedJets += Jet(ptrim, constituents);
      }
      Cut trim_selection = Cuts::abseta < 2.0 && Cuts::pT > 200*GeV && Cuts::massIn(120*GeV, 220*GeV);
      ifilter_select(isortByPt(TrimmedJets), trim_selection);
      if (TrimmedJets.empty())  vetoEvent;


      // SINGLE LEPTON
      bool single_electron=(electrons.size() == 1) && (muons.empty());
      bool single_muon=(muons.size() == 1) && (electrons.empty());

      DressedLepton *lepton = NULL;
      if (single_electron)   lepton = &electrons[0];
      else if (single_muon)  lepton = &muons[0];
      if (!single_electron && !single_muon) vetoEvent;

      //MET
      if (met.pT() < 20*GeV)  vetoEvent;

      //MET+MWT
      double transmass = TransMass(lepton->pt(), lepton->phi(), met.pt(), met.phi());
      if ((met.pT() + transmass) < 60*GeV)   vetoEvent;

      size_t subjets = 0;
      bool btag_hadside=false;
      bool hasHadTopCandidate = false;
      FourMomentum HadTopCandidate;
      for (const Jet& rc_jet : TrimmedJets) {
        FourMomentum rc_jet_mom = rc_jet.mom();
        if (rc_jet_mom.pt() < 350*GeV)   continue;
        double dPhi_lepJet = fabs(deltaPhi(rc_jet_mom.phi(), lepton->phi()));
        if (dPhi_lepJet < 1.)  continue;
        if (rc_jet.pseudojet().user_index()) {
          btag_hadside=true;
        }
        HadTopCandidate = momentum(rc_jet);
        subjets = rc_jet.constituents().size();
        hasHadTopCandidate = true;
        break;
      }
      if (!hasHadTopCandidate)  vetoEvent;

      Jets LepTopCandidates = filter_discard(jets, [&](const Jet& j) {
          return deltaR(j, HadTopCandidate) < 1.5 || deltaR(j, *lepton) > 2.0;
        });
      if (LepTopCandidates.empty()) vetoEvent;

      FourMomentum ltop;
      bool btag_lepside=false;
      for (const Jet& jet : LepTopCandidates) {
        if (jet.bTagged(Cuts::pT > 5*GeV)) {
          btag_lepside = true;
          ltop = jet.mom();
          break;
        }
      }
      if (!btag_hadside && !btag_lepside)  vetoEvent;
      if (!btag_lepside)   ltop = LepTopCandidates[0].momentum();
      double pz = computeneutrinoz(lepton->momentum(), met);
      FourMomentum neutrino( sqrt(sqr(met.px()) + sqr(met.py()) + sqr(pz)), met.px(), met.py(), pz);
      FourMomentum LeptonicTop = lepton->momentum() + neutrino + ltop;
      FourMomentum HadronicTop = HadTopCandidate;
      FourMomentum pttbar = HadronicTop + LeptonicTop;

      Vector3 z_versor(0,0,1);
      Vector3 vpseudotophadron = HadronicTop.vector3();
      Vector3 vpseudotoplepton = LeptonicTop.vector3();
      // Variables
      double ystar = (HadronicTop.pt() > LeptonicTop.pt()) ? 0.5 * (HadronicTop.rap()-LeptonicTop.rap()) : 0.5*(LeptonicTop.rap()-HadronicTop.rap());
      double chi_ttbar = exp(2 * fabs(ystar));
      double pt_leading = (HadronicTop.pt() > LeptonicTop.pt()) ? HadronicTop.pt() : LeptonicTop.pt();
      double pt_subleading = (HadronicTop.pt() > LeptonicTop.pt()) ? LeptonicTop.pt() : HadronicTop.pt();
      double HT_ttbar = HadronicTop.pt() + LeptonicTop.pt();
      double absPout_lep = fabs(vpseudotoplepton.dot((vpseudotophadron.cross(z_versor))/(vpseudotophadron.cross(z_versor).mod())));
      size_t extrajet = smallRjets.size() - subjets -1;

      size_t new_subjets_multi = TransformExtrajetMultiplicity_boosted(subjets);
      size_t new_extrajet_multi = TransformExtrajetMultiplicity_boosted(extrajet);
      size_t new_extrajet_multi_pttop = TransformJetMultiplicity_pttop(extrajet);
      size_t new_extrajet_multi_ptttbar = TransformJetMultiplicity_ptttbar(extrajet);
      size_t new_extrajet_multi_mttbar = TransformJetMultiplicity_mttbar(extrajet);

      _h_multi["boosted_rc_pttop_etatop_multi"].fill(HadronicTop.absrap(), HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_pttop_etattbar_multi"].fill(pttbar.absrap(), HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_pttop_ptttbar_multi"].fill(pttbar.pt()/GeV, HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_pttop_mttbar_multi"].fill(pttbar.mass()/GeV, HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_mttbar_etattbar_multi"].fill(pttbar.absrap(), pttbar.mass()/GeV);
      _h_multi["boosted_rc_mttbar_ptttbar_multi"].fill(pttbar.pt()/GeV, pttbar.mass()/GeV);
      _h_multi["boosted_rc_mttbar_HT_multi"].fill(HT_ttbar, pttbar.mass()/GeV);

      _h_multi["boosted_rc_pttop_extrajet_multi"].fill(new_extrajet_multi_pttop, HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_ptttbar_extrajet_multi"].fill(new_extrajet_multi_ptttbar, pttbar.pt()/GeV);
      _h_multi["boosted_rc_mttbar_extrajet_multi"].fill(new_extrajet_multi_mttbar, pttbar.mass()/GeV);

      _h_multi["boosted_rc_pttop_etatop_multi_norm"].fill(HadronicTop.absrap(), HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_pttop_etattbar_multi_norm"].fill(pttbar.absrap(), HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_pttop_ptttbar_multi_norm"].fill(pttbar.pt()/GeV, HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_pttop_mttbar_multi_norm"].fill(pttbar.mass()/GeV, HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_mttbar_etattbar_multi_norm"].fill(pttbar.absrap(), pttbar.mass()/GeV);
      _h_multi["boosted_rc_mttbar_ptttbar_multi_norm"].fill(pttbar.pt()/GeV, pttbar.mass()/GeV);
      _h_multi["boosted_rc_mttbar_HT_multi_norm"].fill(HT_ttbar/GeV, pttbar.mass()/GeV);
      _h_multi["boosted_rc_pttop_extrajet_multi_norm"].fill(new_extrajet_multi_pttop, HadronicTop.pt()/GeV);
      _h_multi["boosted_rc_ptttbar_extrajet_multi_norm"].fill(new_extrajet_multi_ptttbar, pttbar.pt()/GeV);
      _h_multi["boosted_rc_mttbar_extrajet_multi_norm"].fill(new_extrajet_multi_mttbar, pttbar.mass()/GeV);

      _h["hadTop_boosted_rc_pt"]->fill(HadronicTop.pt()/GeV);
      _h["hadTop_boosted_rc_y"]->fill(HadronicTop.absrap());
      _h["LeadingTop_boosted_rc_pt"]->fill(pt_leading/GeV);
      _h["SubLeadingTop_boosted_rc_pt"]->fill(pt_subleading/GeV);
      _h["boosted_rc_Pout_lep"]->fill(absPout_lep);
      _h["boosted_rc_chi_tt"]->fill(chi_ttbar);
      _h["boosted_rc_HT"]->fill(HT_ttbar/GeV);
      _h["hadTop_boosted_rc_subjets"]->fill(new_subjets_multi);
      _h["boosted_rc_extrajet"]->fill(new_extrajet_multi+1);
      _h["ttbar_boosted_rc_m"]->fill(pttbar.mass()/GeV);
      _h["ttbar_boosted_rc_pt"]->fill(pttbar.pt()/GeV);
      _h["ttbar_boosted_rc_Rapidity"]->fill(pttbar.absrapidity());

      _h["hadTop_boosted_rc_pt_norm"]->fill(HadronicTop.pt()/GeV);
      _h["hadTop_boosted_rc_y_norm"]->fill(HadronicTop.absrap());
      _h["LeadingTop_boosted_rc_pt_norm"]->fill(pt_leading/GeV);
      _h["SubLeadingTop_boosted_rc_pt_norm"]->fill(pt_subleading/GeV);
      _h["boosted_rc_Pout_lep_norm"]->fill(absPout_lep);
      _h["boosted_rc_chi_tt_norm"]->fill(chi_ttbar);
      _h["boosted_rc_HT_norm"]->fill(HT_ttbar/GeV);
      _h["hadTop_boosted_rc_subjets_norm"]->fill(new_subjets_multi);
      _h["boosted_rc_extrajet_norm"]->fill(new_extrajet_multi+1);
      _h["ttbar_boosted_rc_m_norm"]->fill(pttbar.mass()/GeV);
      _h["ttbar_boosted_rc_pt_norm"]->fill(pttbar.pt()/GeV);
      _h["ttbar_boosted_rc_Rapidity_norm"]->fill(pttbar.absrap()) ;
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

    bool _doBoosted, _doResolved;


    double TransMass(double ptLep, double phiLep, double met, double phiMet) {
      return std::sqrt(2.0*ptLep*met*( 1 - std::cos( phiLep-phiMet ) ) );
    }


    double computeneutrinoz(const FourMomentum& lepton, const FourMomentum& met) const {
      //computing z component of neutrino momentum given lepton and met
      double pzneutrino;
      double m_W = 80.399; // in GeV, given in the paper
      double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.px() + lepton.py() * met.py());
      double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
      double b = -2*k*lepton.pz();
      double c = sqr( lepton.E() ) * sqr( met.pT() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
      if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
      else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
        double absquad[2];
        for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
        if (absquad[0] < absquad[1])  pzneutrino = quad[0];
        else                          pzneutrino = quad[1];
      }
      return pzneutrino;
    }


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
        // HepData entry has dummy "Table of Contents", 
        // so need to offset everything by one unit
        { Histo1DPtr tmp; _h_multi[name].add(doubleDiff_bins[i], doubleDiff_bins[i+1], book(tmp, table+1+i, 1, 1)); }
      }
    }


    void book_hist(string name, size_t table) {
      // HepData entry has dummy "Table of Contents", 
      // so need to offset everything by one unit
      book(_h[name], table+3, 1, 1);
      book(_h[name+"_norm"], table+1, 1, 1);
    }

    size_t TransformJetMultiplicity(size_t jet_n) { return jet_n > 7 ? 7 : jet_n; }

    size_t TransformExtrajetMultiplicity(size_t jet_n) { return jet_n > 6 ? 6 : jet_n; }

    size_t TransformExtrajetMultiplicity_boosted(size_t jet_n) { return jet_n > 4 ? 4 : jet_n; }

    size_t TransformJetMultiplicity_for_ttbar_m(size_t jet_n) { return jet_n > 6 ? 6 : jet_n; }

    size_t TransformJetMultiplicity_pttop(size_t jet_n) {
      if (jet_n >= 0 && jet_n < 2)  return 0;
      if (jet_n == 2)               return 2;
      if (jet_n > 2)                return 3;
      return jet_n;
    }

    size_t TransformJetMultiplicity_ptttbar(size_t jet_n) {
      if (jet_n >= 0 && jet_n < 2)  return 0;
      if (jet_n >= 2)               return 2;
      return jet_n;
    }

    size_t TransformJetMultiplicity_mttbar(size_t jet_n) {
      if (jet_n == 0)  return 0;
      if (jet_n == 1)  return 1;
      if (jet_n >= 2)  return 2;
      return jet_n;
    }

    /// @name Objects that are used by the event selection decisions
    ///@{
    map<string, Histo1DPtr> _h;
    map<string, BinnedHistogram> _h_multi;
    ///@}

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1750330);

}
