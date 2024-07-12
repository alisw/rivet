// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {

  /// Colinear Z + Jets in pp at 13 TeV
  class ATLAS_2022_I2077570 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2022_I2077570);

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options
      _mode = 0;
      if ( getOption("LMODE") == "ELEL") _mode = 1;
      if ( getOption("LMODE") == "MUMU") _mode = 2;

      // AntiKt4TruthWZJets prescription
      // Photons
      FinalState all_photons(Cuts::abspid == PID::PHOTON);
      
      // Muons
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true); // true = use muons from prompt tau decays
      DressedLeptons all_dressed_mu(all_photons, bare_mu, 0.1, Cuts::abseta < 2.5, true);
      
      // Electrons
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true); // true = use electrons from prompt tau decays
      DressedLeptons all_dressed_el(all_photons, bare_el, 0.1, Cuts::abseta < 2.5, true);
      
      //Jet forming
      VetoedFinalState vfs(FinalState(Cuts::abseta < 4.5));
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);
            
      FastJets jet(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
      declare(jet, "Jets");


      // Current definition of leptons + jets
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON);
      PromptFinalState muons(Cuts::abspid == PID::MUON);

      // Kinematic cuts for leptons
      const Cut cuts_lep = Cuts::pT > 25*GeV && Cuts::abseta < 2.5;

      DressedLeptons dressed_electrons(photons, electrons, 0.1, cuts_lep);
      declare(dressed_electrons, "DressedElectrons");

      DressedLeptons dressed_muons(photons, muons, 0.1, cuts_lep);
      declare(dressed_muons, "DressedMuons");


      book(_h["ZpT"],        1, 1, 1);
      book(_h["jetpT"],      2, 1, 1);
      book(_h["NJets"],      3, 1, 1);
      book(_h["NJets500"],   4, 1, 1);
      book(_h["minDR"],      5, 1, 1);
      book(_h["rZJ"],        6, 1, 1);
      book(_h["rZJ_coll"],   7, 1, 1);
      book(_h["rZJ_b2b"],    8, 1, 1);
      book(_h["NJets_coll"], 9, 1, 1);
      book(_h["NJets_b2b"], 10, 1, 1);
      book(_h["HT"],        11, 1, 1);
      book(_h["minDR600"],  12, 1, 1);
      book(_h["NJets600"],  13, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // access fiducial electrons and muons
      const Particle *l1 = nullptr, *l2 = nullptr;
      auto muons = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons();
      auto elecs = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons();

      // lepton-jet overlap removal (note: muons are not included in jet finding)
      // Jets eta < 2.5, pT > 30GeV for overlap removal
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 2.5);
      // Remove all jets within dR < 0.2 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets, muons, 0.2);
      idiscardIfAnyDeltaRLess(jets, elecs, 0.2);
      
      // Remove leptons within dR < 0.4 of a jet
      idiscardIfAnyDeltaRLess(muons, jets, 0.4);
      idiscardIfAnyDeltaRLess(elecs, jets, 0.4);

      int nElecs = elecs.size();
      int nMuons = muons.size();

      if (nElecs + nMuons != 2) vetoEvent; // dilepton cut
      if (_mode == 2 && nMuons !=2) vetoEvent;
      if (_mode == 1 && nElecs !=2) vetoEvent;

      if (elecs.size()==2){
        l1=&elecs[0];
        l2=&elecs[1];
      }
      else if (muons.size()==2){
        l1=&muons[0];
        l2=&muons[1];
      }
      else vetoEvent;

      // if not e+e- or mu+mu- pair, veto
      if (l1->pid() + l2->pid() !=0) vetoEvent;

      // Dilepton selection :: Z mass peak.
      FourMomentum ll = (l1->mom() + l2->mom());
      double Zm  = ll.mass(); 
      if ( !inRange(Zm/GeV, 71.0, 111.0) ) vetoEvent;
      double ZpT = ll.pT();
      
      // Calculate the observables
      double jet0pT = 0.;
      double cljetpT = 0.;
      double minDR = 99.;

      // Require jets to be above 100GeV
      ifilter_select(jets, Cuts::pT > 100*GeV);
      double HTjet = sum(jets, Kin::pT, 0.);
      for (const Jet& j : jets) {
        // find minDR and closest jet to Z boson, only with 100GeV+ jets
        double dr = deltaR(j, ll ,RAPIDITY);
        if (dr < minDR) {
          minDR = dr;
          cljetpT = j.pT();
        }
      }
      const size_t Njets = jets.size();

      // NJets >= 1      
      if (Njets < 1) vetoEvent;
      // Exclusive NJets, jet pT > 100 GeV

      _h["NJets"]->fill(Njets);

      // Z pT
      _h["ZpT"]->fill(ZpT/GeV);
      
      //Leading jet pT 
      jet0pT = jets[0].pT();
      _h["jetpT"]->fill(jet0pT/GeV);

      // HT
      double HT = HTjet + l1->pT() + l2->pT();
      _h["HT"]->fill(HT/GeV);

      // HTJet > 600 GeV selection
      if (HTjet >= 600.) {
        double minDR_HTjet = 99.;
        // closest jet to Z
        for (const Jet& j : jets) {
          const double dr = deltaR(j, ll, RAPIDITY);
          if (dr < minDR_HTjet)  minDR_HTjet = dr;
        }
        // Fill histograms of HTjet > 600 GeV
        _h["NJets600"]->fill(Njets);
        _h["minDR600"]->fill(minDR_HTjet);
      } // end of HTjet > 600 GeV

      // Our high pT phase-space
      if (jet0pT/GeV < 500.0) vetoEvent;
      
      // Exclusive NJets, jet pT > 500 GeV
      _h["NJets500"]->fill(Njets);
      
      // Z/J pT ratio
      _h["rZJ"]->fill(ZpT/cljetpT);
      _h["minDR"]->fill(minDR);
      
      // Phase space with DR<1.4       
      if (minDR < 1.4) {
        _h["NJets_coll"]->fill(Njets);
        _h["rZJ_coll"]->fill(ZpT/cljetpT);
      // Phase space with DR>2.0
      } else if (minDR >2.0) {
        _h["NJets_b2b"]->fill(Njets);
        _h["rZJ_b2b"]->fill(ZpT/cljetpT);
      }

    }

    void finalize() {
      double xsec = crossSectionPerEvent()/picobarn;
      // Analysis measures Z->ll(ee or mumu)
      // If file contains both Z->ee and Z->mumu, divide xs by 2
      if (_mode == 0) xsec *= 0.5;
      // Normalize to cross section.
      scale(_h, xsec);
    } // end of finalize

    //@}

    // define histograms
    size_t _mode;
    map<string, Histo1DPtr> _h;
 
  };

  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2077570);
}
