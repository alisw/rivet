#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


class ATLAS_2018_I1705857 : public Analysis {
 public:
   /// Constructor
   /// @brief ttbb at 13 TeV
   DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2018_I1705857);

    void init() {
      // Eta ranges
      Cut eta_full = (Cuts::abseta < 5.0);
      // Lepton cuts
      Cut lep_cuts25 = (Cuts::abseta < 2.5) && (Cuts::pT >= 25*GeV);
      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      PromptFinalState photons(eta_full && Cuts::abspid == PID::PHOTON, true);

      // Projection to find the electrons
      PromptFinalState electrons(eta_full && Cuts::abspid == PID::ELECTRON, true);

      // Projection to find the muons
      PromptFinalState muons(eta_full && Cuts::abspid == PID::MUON, true);

      DressedLeptons dressedelectrons25(photons, electrons, 0.1, lep_cuts25, true);
      DressedLeptons dressedmuons25(photons, muons, 0.1, lep_cuts25, true);

      declare(dressedelectrons25, "elecs");
      declare(dressedmuons25, "muons");

      // From here on we are just setting up the jet clustering
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);

      PromptFinalState jet_photons(eta_full && Cuts::abspid == PID::PHOTON, false);
      DressedLeptons all_dressed_electrons(jet_photons, electrons, 0.1, eta_full, true);
      DressedLeptons all_dressed_muons(jet_photons, muons, 0.1, eta_full, true);

      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(all_dressed_electrons);
      vfs.addVetoOnThisFinalState(all_dressed_muons);
      vfs.addVetoOnThisFinalState(neutrinos);

      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::DECAY_MUONS, JetAlg::DECAY_INVISIBLES);
      declare(jets, "jets");

      // fiducial cross-section histogram
      _histograms["fid_xsec"] = bookHisto1D(1, 1, 1);
      _histograms["fid_xsec_no_ttX"] = bookHisto1D(2, 1, 1);

      _histograms["nbjets_emu"] = bookHisto1D(3, 1, 1);
      _histograms["nbjets_emu_no_ttX"] = bookHisto1D(4, 1, 1);

      // HT
      book_hist("ht_emu", 3);
      book_hist("ht_had_emu", 4);

      book_hist("ht_ljets", 5);
      book_hist("ht_had_ljets", 6);

      // b-jet pTs
      book_hist("lead_bjet_pt_emu",    7);
      book_hist("sublead_bjet_pt_emu", 8);
      book_hist("third_bjet_pt_emu",   9);

      book_hist("lead_bjet_pt_ljets",    10);
      book_hist("sublead_bjet_pt_ljets", 11);
      book_hist("third_bjet_pt_ljets",   12);
      book_hist("fourth_bjet_pt_ljets",  13);

      // leading bb pair
      book_hist("m_bb_leading_emu",  14);
      book_hist("pt_bb_leading_emu", 15);
      book_hist("dR_bb_leading_emu", 16);

      book_hist("m_bb_leading_ljets",  17);
      book_hist("pt_bb_leading_ljets", 18);
      book_hist("dR_bb_leading_ljets", 19);

      // closest bb pair
      book_hist("m_bb_closest_emu",  20);
      book_hist("pt_bb_closest_emu", 21);
      book_hist("dR_bb_closest_emu", 22);

      book_hist("m_bb_closest_ljets",  23);
      book_hist("pt_bb_closest_ljets", 24);
      book_hist("dR_bb_closest_ljets", 25);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      vector<DressedLepton> leptons;
      for (auto &lep : apply<DressedLeptons>(event, "muons").dressedLeptons()) { leptons.push_back(lep); }
      for (auto &lep : apply<DressedLeptons>(event, "elecs").dressedLeptons()) { leptons.push_back(lep); }

      const Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      for (const auto& jet : jets) {
        ifilter_discard(leptons, [&](const DressedLepton& lep) { return deltaR(jet, lep) < 0.4; });
      }

      Jets bjets;
      for (const Jet& jet : jets) {
        if (jet.bTagged(Cuts::pT >= 5*GeV))  bjets += jet;
      }

      size_t njets = jets.size();
      size_t nbjets = bjets.size();

      // Evaluate basic event selection
      bool pass_ljets = (leptons.size() == 1 && leptons[0].pT() > 27*GeV);

      bool pass_emu =
        // 2 leptons > 27 GeV
        (leptons.size() == 2) &&
        (leptons[0].pT() > 27*GeV && leptons[1].pT() > 27*GeV) &&
        // emu events
        ((leptons[0].abspid() == 11 && leptons[1].abspid() == 13) ||
         (leptons[0].abspid() == 13 && leptons[1].abspid() == 11)) &&
        // opposite charge
        (leptons[0].charge() != leptons[1].charge());

      // If we don't have exactly 1 or 2 leptons then veto the event
      if (!pass_emu && !pass_ljets)  vetoEvent;

      if (pass_emu) {
        if (nbjets >= 2)  fill("nbjets_emu", nbjets - 1, weight);
        if (nbjets >= 3)  fill("fid_xsec", 1, weight);
        if (nbjets >= 4)  fill("fid_xsec", 2, weight);
      }

      if (pass_ljets) {
        if (nbjets >= 3 && njets >= 5)  fill("fid_xsec", 3, weight);
        if (nbjets >= 4 && njets >= 6)  fill("fid_xsec", 4, weight);
      }

      if (pass_emu && (nbjets < 3 || njets < 3))    vetoEvent;
      if (pass_ljets && (nbjets < 4 || njets < 6))  vetoEvent;

      double hthad = sum(jets, pT, 0.0);
      double ht = sum(leptons, pT, hthad);

      FourMomentum jsum = bjets[0].momentum() + bjets[1].momentum();
      double dr_leading = deltaR(bjets[0], bjets[1]);

      size_t ind1, ind2; double mindr = 999.;
      for (size_t i = 0; i < bjets.size(); ++i) {
        for (size_t j = 0; j < bjets.size(); ++j) {
          if (i == j)  continue;
          double dr = deltaR(bjets[i], bjets[j]);
          if (dr < mindr) {
            ind1 = i;
            ind2 = j;
            mindr = dr;
          }
        }
      }

      FourMomentum bb_closest = bjets[ind1].momentum() + bjets[ind2].momentum();
      double dr_closest = deltaR(bjets[ind1], bjets[ind2]);

      if (pass_ljets) {
        // b-jet pTs
        fill("lead_bjet_pt_ljets", bjets[0].pT()/GeV, weight);
        fill("sublead_bjet_pt_ljets", bjets[1].pT()/GeV, weight);
        fill("third_bjet_pt_ljets", bjets[2].pT()/GeV, weight);

        if (nbjets >= 4)  fill("fourth_bjet_pt_ljets", bjets[3].pT()/GeV, weight);

        // HT
        fill("ht_ljets", ht/GeV, weight);
        fill("ht_had_ljets", hthad/GeV, weight);

        // leading bb pair
        fill("m_bb_leading_ljets", jsum.mass()/GeV, weight);
        fill("pt_bb_leading_ljets", jsum.pT()/GeV, weight);
        fill("dR_bb_leading_ljets", dr_leading, weight);

        // closest bb pair
        fill("m_bb_closest_ljets", bb_closest.mass()/GeV, weight);
        fill("pt_bb_closest_ljets", bb_closest.pT()/GeV, weight);
        fill("dR_bb_closest_ljets", dr_closest, weight);
      }
      if (pass_emu) {
        // b-jet pTs
        fill("lead_bjet_pt_emu", bjets[0].pT()/GeV, weight);
        fill("sublead_bjet_pt_emu", bjets[1].pT()/GeV, weight);
        fill("third_bjet_pt_emu", bjets[2].pT()/GeV, weight);

        // HT
        fill("ht_emu", ht/GeV, weight);
        fill("ht_had_emu", hthad/GeV, weight);

        // leading bb pair
        fill("m_bb_leading_emu", jsum.mass()/GeV, weight);
        fill("pt_bb_leading_emu", jsum.pT()/GeV, weight);
        fill("dR_bb_leading_emu", dr_leading, weight);

        // closest bb pair
        fill("m_bb_closest_emu", bb_closest.mass()/GeV, weight);
        fill("pt_bb_closest_emu", bb_closest.pT()/GeV, weight);
        fill("dR_bb_closest_emu", dr_closest, weight);
      }
    }

    void finalize() {
      // Normalise all histograms
      const double sf = crossSection() / femtobarn / sumOfWeights();
      for (auto const& h : _histograms) {
        scale(h.second, sf);
        if (h.first.find("fid_xsec") != string::npos)  continue;
        normalize(h.second, 1.0);
      }
    }

    void fill(const string& name, const double value, const double weight) {
      _histograms[name]->fill(value, weight);
      _histograms[name + "_no_ttX"]->fill(value, weight);
    }

    void book_hist(const std::string& name, unsigned int d) {
      // ttX substracted histograms are even numbers
      _histograms[name] = bookHisto1D((d * 2) - 1, 1, 1);
      _histograms[name + "_no_ttX"] = bookHisto1D(d * 2, 1, 1);
    }

    private:
      map<std::string, Histo1DPtr> _histograms;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2018_I1705857);
}
