// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {


  class ATLAS_2019_I1746286 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1746286);

    /// Book histograms and initialise projections before the run
    void init() {
      // Set up projections
      const FinalState fs(Cuts::abseta < 4.5);

      /// Get electrons from truth record
      FinalState elec_fs(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(elec_fs, "ELEC_FS");

      /// Get muons which pass the initial kinematic cuts:
      FinalState muon_fs(Cuts::abspid == PID::MUON && Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      declare(muon_fs, "MUON_FS");

      // get b-hadrons
      declare(HeavyHadrons(Cuts::pT > 5*GeV), "BHadrons");

      UnstableParticles k0_fs(Cuts::abspid == PID::K0S && Cuts::abseta < 2.5 && Cuts::E > 1*GeV);
      declare(k0_fs, "K0_FS");

      UnstableParticles lambda_fs(Cuts::abspid == PID::LAMBDA && Cuts::abseta < 2.5 && Cuts::E > 1*GeV);
      declare(lambda_fs, "LAMBDA_FS");

      // Final state used as input for jet-finding.
      // We include everything except the muons and neutrinos
      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "JETS");

      // Book histograms
      book(_h["b_k0_pt"],  1, 1, 1);
      book(_h["b_k0_x"],   2, 1, 1);
      book(_h["b_k0_e"],   3, 1, 1);
      book(_h["b_k0_eta"], 4, 1, 1);
      book(_h["b_k0_n"],   5, 1, 1);

      book(_h["j_k0_pt"],  6, 1, 1);
      book(_h["j_k0_x"],   7, 1, 1);
      book(_h["j_k0_e"],   8, 1, 1);
      book(_h["j_k0_eta"], 9, 1, 1);
      book(_h["j_k0_n"],  10, 1, 1);

      book(_h["out_k0_pt"],  11, 1, 1);
      book(_h["out_k0_e"],   12, 1, 1);
      book(_h["out_k0_eta"], 13, 1, 1);
      book(_h["out_k0_n"],   14, 1, 1);

      book(_h["all_k0_pt"],  15, 1, 1);
      book(_h["all_k0_e"],   16, 1, 1);
      book(_h["all_k0_eta"], 17, 1, 1);
      book(_h["all_k0_n"],   18, 1, 1);

      book(_h["all_l_pt"],  19, 1, 1);
      book(_h["all_l_e"],   20, 1, 1);
      book(_h["all_l_eta"], 21, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// Get the various sets of final state particles
      const Particles& elecFS = apply<FinalState>(event, "ELEC_FS").particlesByPt();
      const Particles& muonFS = apply<FinalState>(event, "MUON_FS").particlesByPt();
      const Particles& k0FS = apply<UnstableParticles>(event, "K0_FS").particlesByPt();
      const Particles& lambdaFS = apply<UnstableParticles>(event, "LAMBDA_FS").particlesByPt();

      // Get all jets with pT > 7 GeV (ATLAS standard jet collection)
      Jets jets = apply<FastJets>(event, "JETS").jetsByPt(7*GeV);

      // Keep any jets that pass the pt cut
      Jets good_jets = filter_select(jets, Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      // Remove jets too close to an electron
      idiscardIfAnyDeltaRLess(good_jets, elecFS, 0.2);

      // Classify the event type
      const size_t nElec = elecFS.size();
      const size_t nMuon = muonFS.size();
      bool isDilepton = false;
      if (nElec == 2 && nMuon == 0) {
        if (charge(elecFS[0]) != charge(elecFS[1])) isDilepton = true;
      } else if (nElec == 1 && nMuon == 1) {
        if (charge(elecFS[0]) != charge(muonFS[0])) isDilepton = true;
      } else if (nElec == 0 && nMuon == 2) {
        if (charge(muonFS[0]) != charge(muonFS[1])) isDilepton = true;
      }
      const bool isGoodEvent = (isDilepton && good_jets.size() >= 2);
      if (!isGoodEvent) vetoEvent;

      // Select b-hadrons
      const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();

      // Select b-jets as those containing a b-hadron
      Jets bjets = discardIfAnyDeltaRLess(good_jets, bHadrons, 0.3);

      size_t n_k0_all = 0;
      size_t n_k0_out = 0;
      size_t n_k0_b = 0;
      size_t n_k0_j = 0;

      size_t n_k0_all_visible = 0;
      size_t n_k0_out_visible = 0;
      size_t n_k0_b_visible = 0;
      size_t n_k0_j_visible = 0;

      bool isVisible = false;
      // Loop over all K0s particles
      for (const Particle& k : k0FS) {
        if (k.hasStableDescendantWith(Cuts::pid == PID::PIPLUS)) isVisible = true;

        n_k0_all += 1;
        if (isVisible) n_k0_all_visible += 1;
        _h["all_k0_pt"]->fill(k.pT()/GeV);
        _h["all_k0_eta"]->fill(k.abseta());
        _h["all_k0_e"]->fill(k.E()/GeV);

        bool isJetAssoc = false, isBjet = false;
        double minDeltaR = 1000., jetAssocE = 0.;

        for (const Jet& j : good_jets) {
          const double k0_jetdR = deltaR(j, k);
          if (k0_jetdR < 0.4 && k0_jetdR < minDeltaR) { 
            isJetAssoc = true;
            minDeltaR = k0_jetdR;
            jetAssocE = j.E();
            isBjet = any(bHadrons, DeltaRLess(j, 0.3));
          }
        }

        // K0s not associated to jets
        if (!isJetAssoc){
          n_k0_out += 1;
          if(isVisible) n_k0_out_visible += 1;
          _h["out_k0_pt"]->fill(k.pT()/GeV);
          _h["out_k0_eta"]->fill(k.abseta());
          _h["out_k0_e"]->fill(k.E()/GeV);
        }

        //K0s associated to b-jets
        if(isJetAssoc && isBjet){
          n_k0_b += 1;
          if(isVisible) n_k0_b_visible += 1;
          _h["b_k0_pt"]->fill(k.pT()/GeV);
          _h["b_k0_eta"]->fill(k.abseta());
          _h["b_k0_e"]->fill(k.E()/GeV);
          _h["b_k0_x"]->fill(k.E()/jetAssocE);
        }

        //K0s associated to non b-jets
        if(isJetAssoc && !isBjet){
          n_k0_j += 1;
          if(isVisible) n_k0_j_visible += 1;
          _h["j_k0_pt"]->fill(k.pT()/GeV);
          _h["j_k0_eta"]->fill(k.abseta());
          _h["j_k0_e"]->fill(k.E()/GeV);
          _h["j_k0_x"]->fill(k.E()/jetAssocE);
        }
      }


      // K0s multiplicities
      _h["all_k0_n"]->fill(n_k0_all_visible);
      _h["out_k0_n"]->fill(n_k0_out_visible);
      _h["b_k0_n"]->fill(n_k0_b_visible);
      _h["j_k0_n"]->fill(n_k0_j_visible);
    

      // Loop over all Lambda particles
      //size_t n_lambda_all = 0;
      for(const Particle& l : lambdaFS) {
        //n_lambda_all += 1;
        _h["all_l_pt"]->fill(l.pT()/GeV);
        _h["all_l_eta"]->fill(l.abseta());
        _h["all_l_e"]->fill(l.E()/GeV);
      }
    }


    // Histogram normalization to the number of events passing the cuts
    void finalize(){
      const double sf = 1.0 / _h["all_k0_n"]->sumW();
      for (auto& hist : _h) { scale(hist.second, sf); }
    }

  private:

    // Counters
    map<string, Histo1DPtr> _h;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1746286);
}
