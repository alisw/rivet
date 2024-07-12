#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {

  /// ttbar dilepton differential cross-sections in pp collisions at 13 TeV
  class CMS_2018_I1703993 : public Analysis {  //
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1703993);  //

    void init() {
      // Parton level top quarks dilepton e/mu channels only
      const bool acceptTauDecays = false;
      declare(PartonicTops(PartonicTops::DecayMode::E_MU, acceptTauDecays),
              "PartonTops");  // Partonic top decaying to e or mu

      // Build particle level tops starting from FinalState
      const FinalState fs(Cuts::pT > 0. && Cuts::abseta < 6.);

      // Neutrinos
      InvisibleFinalState neutrinos(true, true, true);
      declare(neutrinos, "Neutrinos");

      // Projection for electrons and muons
      Cut lepton_cut = Cuts::pt > 20 * GeV && Cuts::abseta < 2.4;

      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      PromptFinalState prompt_leptons(charged_leptons, false);
      PromptFinalState prompt_photons(photons, true);

      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, lepton_cut, /*useDecayPhotons*/ true);
      declare(dressed_leptons, "DressedLeptons");

      // Projection for jets
      VetoedFinalState fs_jets(fs);
      fs_jets.addVetoOnThisFinalState(dressed_leptons);
      fs_jets.vetoNeutrinos();
      declare(FastJets(fs_jets, FastJets::ANTIKT, 0.4), "ak4jets");

      // Book hists for particle level, normalized values
      book(_h["top_pt_norm"], "d07-x01-y01");       //  table 15, normalized values
      book(_h["atop_pt_norm"], "d15-x01-y01");      //  table 16, normalized values
      book(_h["top1_pt_norm"], "d23-x01-y01");      //  table 17, normalized values
      book(_h["top2_pt_norm"], "d31-x01-y01");      //  table 18, normalized values
      book(_h["top_pt_ttrf_norm"], "d39-x01-y01");  //  table 19, normalized values
      book(_h["top_y_norm"], "d47-x01-y01");        //  table 20, normalized values
      book(_h["atop_y_norm"], "d55-x01-y01");       //  table 21, normalized values
      book(_h["top1_y_norm"], "d63-x01-y01");       //  table 22, normalized values
      book(_h["top2_y_norm"], "d71-x01-y01");       //  table 23, normalized values
      book(_h["tt_pt_norm"], "d79-x01-y01");        //  table 24, normalized values
      book(_h["tt_y_norm"], "d87-x01-y01");         //  table 25, normalized values
      book(_h["tt_m_norm"], "d95-x01-y01");         //  table 26, normalized values
      book(_h["tt_dy_norm"], "d103-x01-y01");       //  table 27, normalized values
      book(_h["tt_dphi_norm"], "d111-x01-y01");     //  table 28, normalized values
      book(_h["lep_pt_norm"], "d115-x01-y01");      //  table 29, normalized values,lepton
      book(_h["alep_pt_norm"], "d119-x01-y01");     //  table 30, normalized values,antilepton
      book(_h["lep1_pt_norm"], "d123-x01-y01");     //  table 31, normalized values,leading
      book(_h["lep2_pt_norm"], "d127-x01-y01");     //  table 32, normalized values,trailing
      book(_h["lep_eta_norm"], "d131-x01-y01");     //  table 33, normalized values, eta lepton
      book(_h["alep_eta_norm"], "d135-x01-y01");    //  table 34, normalized values, eta antilepton
      book(_h["lep1_eta_norm"], "d139-x01-y01");    //  table 35, normalized values, eta l leading
      book(_h["lep2_eta_norm"], "d143-x01-y01");    //  table 36, normalized values, eta l trailing
      book(_h["ll_pt_norm"], "d147-x01-y01");       // table 37, normalized values
      book(_h["ll_m_norm"], "d151-x01-y01");        // table 38, normalized values
      book(_h["ll_dphi_norm"], "d155-x01-y01");     // table 39, normalized values
      book(_h["ll_deta_norm"], "d159-x01-y01");     // table 40, normalized values
      book(_h["jet_sz_norm"], "d187-x01-y01");      //  table 41(47), normalized values
      book(_h["jet1_pt_norm"], "d163-x01-y01");     //  table 42, normalized values, leading
      book(_h["jet2_pt_norm"], "d167-x01-y01");     //  table 43, normalized values
      book(_h["b1_eta_norm"], "d171-x01-y01");      //  table 44(43), normalized values, leading
      book(_h["b2_eta_norm"], "d175-x01-y01");      //  table 45(44), normalized values, trailing
      book(_h["bb_pt_norm"], "d179-x01-y01");       //  table 46(45), normalized values
      book(_h["bb_m_norm"], "d183-x01-y01");        //  table 47(46), normalized values

      // book hists for particle level, absolute values
      book(_h_abs["top_pt"], "d05-x01-y01");       //  table 15
      book(_h_abs["atop_pt"], "d13-x01-y01");      //  table 16
      book(_h_abs["top1_pt"], "d21-x01-y01");      //  table 17
      book(_h_abs["top2_pt"], "d29-x01-y01");      //  table 18
      book(_h_abs["top_pt_ttrf"], "d37-x01-y01");  //  table 19
      book(_h_abs["top_y"], "d45-x01-y01");        //  table 20
      book(_h_abs["atop_y"], "d53-x01-y01");       //  table 21
      book(_h_abs["top1_y"], "d61-x01-y01");       //  table 22
      book(_h_abs["top2_y"], "d69-x01-y01");       //  table 23
      book(_h_abs["tt_pt"], "d77-x01-y01");        //  table 24
      book(_h_abs["tt_y"], "d85-x01-y01");         //  table 25
      book(_h_abs["tt_m"], "d93-x01-y01");         //  table 26
      book(_h_abs["tt_dy"], "d101-x01-y01");       //  table 27
      book(_h_abs["tt_dphi"], "d109-x01-y01");     //  table 28
      book(_h_abs["lep_pt"], "d113-x01-y01");      //  table 29, lepton
      book(_h_abs["alep_pt"], "d117-x01-y01");     //  table 30, antilepton
      book(_h_abs["lep1_pt"], "d121-x01-y01");     //  table 31, leading
      book(_h_abs["lep2_pt"], "d125-x01-y01");     //  table 32, trailing
      book(_h_abs["lep_eta"], "d129-x01-y01");     //  table 33, eta lepton
      book(_h_abs["alep_eta"], "d133-x01-y01");    //  table 34, eta antilepton
      book(_h_abs["lep1_eta"], "d137-x01-y01");    //  table 35, eta l leading
      book(_h_abs["lep2_eta"], "d141-x01-y01");    //  table 36, eta l trailing
      book(_h_abs["ll_pt"], "d145-x01-y01");       // table 37
      book(_h_abs["ll_m"], "d149-x01-y01");        // table 38
      book(_h_abs["ll_dphi"], "d153-x01-y01");     // table 39
      book(_h_abs["ll_deta"], "d157-x01-y01");     // table 40
      book(_h_abs["jet_sz"], "d185-x01-y01");      //  table 41(47)
      book(_h_abs["jet1_pt"], "d161-x01-y01");     //  table 42(41), leading
      book(_h_abs["jet2_pt"], "d165-x01-y01");     //  table 43(42)
      book(_h_abs["b1_eta"], "d169-x01-y01");      //  table 44(43), leading
      book(_h_abs["b2_eta"], "d173-x01-y01");      //  table 45(44), trailing
      book(_h_abs["bb_pt"], "d177-x01-y01");       //  table 46(45)
      book(_h_abs["bb_m"], "d181-x01-y01");        //  table 47(46)

      //book hists for parton level, normalized values
      book(_h_part["top_pt_norm"], "d03-x01-y01");       //  table 1
      book(_h_part["atop_pt_norm"], "d11-x01-y01");      //  table 2
      book(_h_part["top1_pt_norm"], "d19-x01-y01");      //  table 3
      book(_h_part["top2_pt_norm"], "d27-x01-y01");      //  table 4
      book(_h_part["top_pt_ttrf_norm"], "d35-x01-y01");  //  table 5
      book(_h_part["top_y_norm"], "d43-x01-y01");        //  table 6
      book(_h_part["atop_y_norm"], "d51-x01-y01");       //  table 7
      book(_h_part["top1_y_norm"], "d59-x01-y01");       //  table 8
      book(_h_part["top2_y_norm"], "d67-x01-y01");       //  table 9
      book(_h_part["tt_pt_norm"], "d75-x01-y01");        //  table 10
      book(_h_part["tt_y_norm"], "d83-x01-y01");         //  table 11
      book(_h_part["tt_m_norm"], "d91-x01-y01");         //  table 12
      book(_h_part["tt_dy_norm"], "d99-x01-y01");        //  table 13
      book(_h_part["tt_dphi_norm"], "d107-x01-y01");     //  table 14

      //book hists for parton level, absolute values
      book(_h_part_abs["top_pt"], "d01-x01-y01");       //  table 1
      book(_h_part_abs["atop_pt"], "d09-x01-y01");      //  table 2
      book(_h_part_abs["top1_pt"], "d17-x01-y01");      //  table 3
      book(_h_part_abs["top2_pt"], "d25-x01-y01");      //  table 4
      book(_h_part_abs["top_pt_ttrf"], "d33-x01-y01");  //  table 5
      book(_h_part_abs["top_y"], "d41-x01-y01");        //  table 6
      book(_h_part_abs["atop_y"], "d49-x01-y01");       //  table 7
      book(_h_part_abs["top1_y"], "d57-x01-y01");       //  table 8
      book(_h_part_abs["top2_y"], "d65-x01-y01");       //  table 9
      book(_h_part_abs["tt_pt"], "d73-x01-y01");        //  table 10
      book(_h_part_abs["tt_y"], "d81-x01-y01");         //  table 11
      book(_h_part_abs["tt_m"], "d89-x01-y01");         //  table 12
      book(_h_part_abs["tt_dy"], "d97-x01-y01");        //  table 13
      book(_h_part_abs["tt_dphi"], "d105-x01-y01");     //  table 14
    }

    void analyze(const Event& event) {
      // Parton-level analysis
      const Particles partonTops = apply<ParticleFinder>(event, "PartonTops").particlesByPt();
      if (partonTops.size() == 2) {
        Particle top1 = partonTops[0];
        Particle top2 = partonTops[1];

        Particle top = top1;
        Particle atop = top2;
        if (top.pid() < 0) {
          std::swap(top, atop);
        }

        const FourMomentum tt = top1.momentum() + top2.momentum();
        const LorentzTransform rf = LorentzTransform::mkFrameTransform(tt);
        const FourMomentum top_rf = rf.transform(top);

        //Parton level normalized values
        _h_part["top_pt_norm"]->fill(top.pt()/GeV);
        _h_part["atop_pt_norm"]->fill(atop.pt()/GeV);
        _h_part["top1_pt_norm"]->fill(top1.pt()/GeV);
        _h_part["top2_pt_norm"]->fill(top2.pt()/GeV);
        _h_part["top_pt_ttrf_norm"]->fill(top_rf.pt()/GeV);
        _h_part["top_y_norm"]->fill(top.rapidity());
        _h_part["atop_y_norm"]->fill(atop.rapidity());
        _h_part["top1_y_norm"]->fill(top1.rapidity());
        _h_part["top2_y_norm"]->fill(top2.rapidity());
        _h_part["tt_pt_norm"]->fill(tt.pt()/GeV);
        _h_part["tt_y_norm"]->fill(tt.rapidity());
        _h_part["tt_m_norm"]->fill(tt.mass()/GeV);
        _h_part["tt_dy_norm"]->fill(deltaRap(top.absrap(), atop.absrap(), true));
        _h_part["tt_dphi_norm"]->fill(deltaPhi(top.phi(), atop.phi()));

        //Parton level absolute values
        _h_part_abs["top_pt"]->fill(top.pt()/GeV);
        _h_part_abs["atop_pt"]->fill(atop.pt()/GeV);
        _h_part_abs["top1_pt"]->fill(top1.pt()/GeV);
        _h_part_abs["top2_pt"]->fill(top2.pt()/GeV);
        _h_part_abs["top_pt_ttrf"]->fill(top_rf.pt()/GeV);
        _h_part_abs["top_y"]->fill(top.rapidity());
        _h_part_abs["atop_y"]->fill(atop.rapidity());
        _h_part_abs["top1_y"]->fill(top1.rapidity());
        _h_part_abs["top2_y"]->fill(top2.rapidity());
        _h_part_abs["tt_pt"]->fill(tt.pt()/GeV);
        _h_part_abs["tt_y"]->fill(tt.rapidity());
        _h_part_abs["tt_m"]->fill(tt.mass()/GeV);
        _h_part_abs["tt_dy"]->fill(deltaRap(top.absrap(), atop.absrap(), true));
        _h_part_abs["tt_dphi"]->fill(deltaPhi(top.phi(), atop.phi()));
      }

      //Particle-level analysis
      // Select leptons
      const vector<DressedLepton>& dressedLeptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      if (dressedLeptons.size() != 2)
        vetoEvent;
      sortByPt(dressedLeptons);

      const FourMomentum& lepton1 = dressedLeptons[0].momentum();
      const FourMomentum& lepton2 = dressedLeptons[1].momentum();
      if ((lepton1 + lepton2).mass() < 20*GeV)
        vetoEvent;

      // Select neutrinos
      const Particles neutrinos = apply<InvisibleFinalState>(event, "Neutrinos").particlesByPt();
      if (neutrinos.size() < 2)
        vetoEvent;

      // Select bjets
      const FastJets& fjJets = apply<FastJets>(event, "ak4jets");
      const Jets jets = fjJets.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30 * GeV);
      const Jets bJets = filter_select(jets, hasBTag());

      // There should at least two b jets.
      if (bJets.size() < 2)
        vetoEvent;

      // Construct particle level top
      FourMomentum nu1 = neutrinos[0].momentum();
      FourMomentum nu2 = neutrinos[1].momentum();
      if (std::abs((lepton1 + nu1).mass() - 80.4*GeV) + std::abs((lepton2 + nu2).mass() - 80.4*GeV) >
          std::abs((lepton1 + nu2).mass() - 80.4*GeV) + std::abs((lepton2 + nu1).mass() - 80.4*GeV)) {
        std::swap(nu1, nu2);
      }
      const FourMomentum w1 = lepton1 + nu1;
      const FourMomentum w2 = lepton2 + nu2;

      double dm = 1e9;  // Reset once again for top combination.
      int selB1 = -1, selB2 = -1;
      for (unsigned int i = 0; i < bJets.size(); ++i) {
        const auto& bjet1 = bJets.at(i);
        if (deltaR(bjet1, lepton1) < 0.4)
          continue;
        if (deltaR(bjet1, lepton2) < 0.4)
          continue;

        const double dm1 = std::abs((w1 + bjet1).mass() - 172.5*GeV);
        for (unsigned int j = 0; j < bJets.size(); ++j) {
          if (i == j)
            continue;
          const auto& bjet2 = bJets.at(j);
          if (deltaR(bjet2, lepton1) < 0.4)
            continue;
          if (deltaR(bjet2, lepton2) < 0.4)
            continue;
          const double dm2 = std::abs((w2 + bjet2).mass() - 172.5*GeV);
          const double newDm = dm1 + dm2;

          if (newDm < dm) {
            dm = newDm;
            selB1 = i;
            selB2 = j;
          }
        }
      }

      if (dm >= 1e9)
        vetoEvent;

      FourMomentum bjet1 = bJets[selB1].momentum();
      FourMomentum bjet2 = bJets[selB2].momentum();

      const FourMomentum t1 = w1 + bjet1;
      const FourMomentum t2 = w2 + bjet2;
      const FourMomentum tt = t1 + t2;

      int q1 = dressedLeptons[0].charge();
      int q2 = dressedLeptons[1].charge();
      if (q1 * q2 > 0)
        vetoEvent;

      int idx_lep = (q1 == -1) ? 0 : 1;
      int idx_alep = (q2 == 1) ? 1 : 0;

      FourMomentum top = t1;
      FourMomentum atop = t2;
      if (q1 == -1) {
        std::swap(top, atop);
      }

      FourMomentum top1 = t1;
      FourMomentum top2 = t2;
      if (top1.pt() < top2.pt()) {
        std::swap(top1, top2);
      }

      if (bjet1.pt() < bjet2.pt()) {
        std::swap(bjet1, bjet2);
      }

      const LorentzTransform rf = LorentzTransform::mkFrameTransform(tt);
      const FourMomentum top_rf = rf.transform(top);

      const FourMomentum ll = lepton1 + lepton2;

      const FourMomentum bb = bjet1 + bjet2;

      //particle level normalized values
      _h["top_pt_norm"]->fill(top.pt()/GeV);
      _h["atop_pt_norm"]->fill(atop.pt()/GeV);
      _h["top1_pt_norm"]->fill(top1.pt()/GeV);
      _h["top2_pt_norm"]->fill(top2.pt()/GeV);
      _h["top_pt_ttrf_norm"]->fill(top_rf.pt()/GeV);
      _h["top_y_norm"]->fill(top.rapidity());
      _h["atop_y_norm"]->fill(atop.rapidity());
      _h["top1_y_norm"]->fill(top1.rapidity());
      _h["top2_y_norm"]->fill(top2.rapidity());
      _h["tt_pt_norm"]->fill(tt.pt()/GeV);
      _h["tt_dy_norm"]->fill(deltaRap(top.absrap(), atop.absrap(), true));
      _h["tt_y_norm"]->fill(tt.rapidity());
      _h["tt_m_norm"]->fill(tt.mass()/GeV);
      _h["tt_dphi_norm"]->fill(deltaPhi(top.phi(), atop.phi()));
      _h["lep_pt_norm"]->fill(dressedLeptons[idx_lep].pt()/GeV);
      _h["alep_pt_norm"]->fill(dressedLeptons[idx_alep].pt()/GeV);
      _h["lep1_pt_norm"]->fill(lepton1.pt()/GeV);
      _h["lep2_pt_norm"]->fill(lepton2.pt()/GeV);
      _h["lep_eta_norm"]->fill(dressedLeptons[idx_lep].eta());
      _h["alep_eta_norm"]->fill(dressedLeptons[idx_alep].eta());
      _h["lep1_eta_norm"]->fill(lepton1.eta());
      _h["lep2_eta_norm"]->fill(lepton2.eta());
      _h["ll_pt_norm"]->fill(ll.pt()/GeV);
      _h["ll_m_norm"]->fill(ll.mass()/GeV);
      _h["ll_dphi_norm"]->fill(deltaPhi(dressedLeptons[idx_lep], dressedLeptons[idx_alep]));
      _h["ll_deta_norm"]->fill(deltaEta(dressedLeptons[idx_lep].abseta(), dressedLeptons[idx_alep].abseta(), true));
      _h["jet_sz_norm"]->fill(jets.size());
      _h["jet1_pt_norm"]->fill(bjet1.pt()/GeV);
      _h["jet2_pt_norm"]->fill(bjet2.pt()/GeV);
      _h["b1_eta_norm"]->fill(bjet1.eta());
      _h["b2_eta_norm"]->fill(bjet2.eta());
      _h["bb_pt_norm"]->fill(bb.pt()/GeV);
      _h["bb_m_norm"]->fill(bb.mass()/GeV);

      // absolute values
      _h_abs["top_pt"]->fill(top.pt()/GeV);
      _h_abs["atop_pt"]->fill(atop.pt()/GeV);
      _h_abs["top1_pt"]->fill(top1.pt()/GeV);
      _h_abs["top2_pt"]->fill(top2.pt()/GeV);
      _h_abs["top_pt_ttrf"]->fill(top_rf.pt()/GeV);
      _h_abs["top_y"]->fill(top.rapidity());
      _h_abs["atop_y"]->fill(atop.rapidity());
      _h_abs["top1_y"]->fill(top1.rapidity());
      _h_abs["top2_y"]->fill(top2.rapidity());
      _h_abs["tt_pt"]->fill(tt.pt()/GeV);
      _h_abs["tt_y"]->fill(tt.rapidity());
      _h_abs["tt_m"]->fill(tt.mass()/GeV);
      _h_abs["tt_dy"]->fill(deltaRap(top.absrap(), atop.absrap(), true));
      _h_abs["tt_dphi"]->fill(deltaPhi(top.phi(), atop.phi()));
      _h_abs["lep_pt"]->fill(dressedLeptons[idx_lep].pt()/GeV);
      _h_abs["alep_pt"]->fill(dressedLeptons[idx_alep].pt()/GeV);
      _h_abs["lep1_pt"]->fill(lepton1.pt()/GeV);
      _h_abs["lep2_pt"]->fill(lepton2.pt()/GeV);
      _h_abs["lep_eta"]->fill(dressedLeptons[idx_lep].eta());
      _h_abs["alep_eta"]->fill(dressedLeptons[idx_alep].eta());
      _h_abs["lep1_eta"]->fill(lepton1.eta());
      _h_abs["lep2_eta"]->fill(lepton2.eta());
      _h_abs["ll_pt"]->fill(ll.pt()/GeV);
      _h_abs["ll_m"]->fill(ll.mass()/GeV);
      _h_abs["ll_dphi"]->fill(deltaPhi(dressedLeptons[idx_lep], dressedLeptons[idx_alep]));
      _h_abs["ll_deta"]->fill(deltaEta(dressedLeptons[idx_lep].abseta(), dressedLeptons[idx_alep].abseta(), true));
      _h_abs["jet_sz"]->fill(jets.size());
      _h_abs["jet1_pt"]->fill(bjet1.pt()/GeV);
      _h_abs["jet2_pt"]->fill(bjet2.pt()/GeV);
      _h_abs["b1_eta"]->fill(bjet1.eta());
      _h_abs["b2_eta"]->fill(bjet2.eta());
      _h_abs["bb_pt"]->fill(bb.pt()/GeV);
      _h_abs["bb_m"]->fill(bb.mass()/GeV);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
      scale(_h_abs, crossSection() / picobarn / sumOfWeights());

      normalize(_h_part);
      scale(_h_part_abs,
            crossSection() / picobarn / sumOfWeights() /
                0.04553956);  // BR correction for E_MU: 0.04553956  //BR correction for MUON: 0.01129969
    }

  private:
    map<string, Histo1DPtr> _h;
    map<string, Histo1DPtr> _h_abs;
    map<string, Histo1DPtr> _h_part;
    map<string, Histo1DPtr> _h_part_abs;
  };

  RIVET_DECLARE_PLUGIN(CMS_2018_I1703993);

}  // namespace Rivet
