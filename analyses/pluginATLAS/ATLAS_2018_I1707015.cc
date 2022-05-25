// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  /// @brief: ttbar + gamma at 13 TeV
  class ATLAS_2018_I1707015 : public Analysis {
  public:


    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1707015);


    // Book histograms and initialise projections before the run
    void init() {


      // Set default running mode to 3 (all)
      _mode = 3;

      // Get running mode
      if ( getOption("LMODE") == "SINGLE" )   _mode = 1;
      if ( getOption("LMODE") == "DILEPTON" ) _mode = 2;
      if ( getOption("LMODE") == "ALL" )      _mode = 3;

      // All final state particles
      const FinalState fs;

      // Charged particles for signal photon isolation
      ChargedFinalState cfs(fs);
      declare(cfs, "CFS");

      // Signal photons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON && Cuts::pT > 20*GeV && Cuts::abseta < 2.37, true);
      declare(photons, "Photons");

      // Leptons
      PromptFinalState leptons(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON, true);

      // Dress the leptons
      FinalState dressPhotons(Cuts::abspid == PID::PHOTON);
      Cut lepCuts = (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV);
      DressedLeptons dressedLeptons(dressPhotons, leptons, 0.1, lepCuts);
      declare(dressedLeptons, "Leptons");

      // Jet alg input
      VetoedFinalState vfs(fs);

      // Remove prompt invisibles from jet input
      VetoedFinalState invis_fs(fs);
      invis_fs.addVetoOnThisFinalState(VisibleFinalState(fs));
      PromptFinalState invis_pfs = PromptFinalState(invis_fs, true);
      vfs.addVetoOnThisFinalState(invis_pfs);

      // Remove prompt dressed muons (muons + associated photons) from jet input
      PromptFinalState muons(Cuts::abspid == PID::MUON, true);
      DressedLeptons dressedmuons(dressPhotons, muons, 0.1);
      vfs.addVetoOnThisFinalState(dressedmuons);

      // Jet clustering
      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(jets, "Jets");

      // Book histograms
      if ( _mode == 1 or _mode == 3 ){  // Single lepton channel
        book(_h["sl_ph_pt"],   3,1,1);  // photon pT
        book(_h["sl_ph_eta"],  4,1,1);  // photon eta
        book(_h["sl_ph_l_dR"], 5,1,1);  // photon-lepton dR
      }
      if ( _mode == 2 or _mode == 3 ){  // Dilepton channel
        book(_h["dl_ph_pt"],   6,1,1);  // photon pT
        book(_h["dl_ph_eta"],  7,1,1);  // photon eta
        book(_h["dl_ph_l_dR"], 8,1,1);  // min photon-lepton dR
        book(_h["dl_l_dEta"],  9,1,1);  // lepton-lepton dEta
        book(_h["dl_l_dPhi"],  10,1,1); // lepton-lepton dPhi
      }
    }


    void analyze(const Event& event) {

      // Fetch objects
      const vector<DressedLepton>& leptons = apply<DressedLeptons>(event, "Leptons").dressedLeptons();
      Particles photons = applyProjection<PromptFinalState>(event, "Photons").particles();
      ChargedFinalState charged = apply<ChargedFinalState>(event, "CFS");
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);

      // Immediate veto on events without one good photon
      if ( photons.size() != 1 ) vetoEvent;
      const Particle& photon = photons[0];

      // Veto event if photon too close to a lepton
      for (const DressedLepton& lep : leptons) {
        if ( deltaR(lep, photon) < 1.0 ) vetoEvent;
      }

      // Overlap removel of jets near leptons
      for (const DressedLepton& lep : leptons) {
        ifilter_discard(jets, deltaRLess(lep, 0.4));
      }

      // Overlap removel of jets near isolated photon
      double conePt = 0.0;
      Particles photSurround = charged.particles(deltaRLess(photon, 0.3));
      for (const Particle& p : photSurround) {
        conePt += p.pt();
      }
      if ( conePt / photon.pT() < 0.1 ) ifilter_discard(jets, deltaRLess(photon, 0.4) );

      // Veto event if photon too close to good jets
      for (const Jet& jet : jets) {
        if ( deltaR(jet, photon) < 0.4 ) vetoEvent;
      }

      // Require at least one bjet
      unsigned int nbjets = 0;
      for (const Jet& jet : jets) {
        if ( jet.bTagged(Cuts::pT > 5*GeV) ) nbjets += 1;
      }
      if ( nbjets == 0 ) vetoEvent;

      // Fill histos in single lepton channel
      if ( _mode == 1 || _mode == 3 ) {
        if ( jets.size() >= 4 && leptons.size() == 1 ) {
          _h["sl_ph_pt"]->fill(   photon.pT()/GeV );
          _h["sl_ph_eta"]->fill(  photon.abseta() );
          _h["sl_ph_l_dR"]->fill( deltaR(leptons[0], photon) );
        }
      }

      // Fill histos in dilepton channel
      if ( _mode == 2 || _mode == 3 ) {
        if ( jets.size() >= 2 && leptons.size() == 2 ) {
          double deltaRNew;
          double deltaRMin = 999.0;
          for (const DressedLepton& lep : leptons) {
            deltaRNew = deltaR(lep, photon);
            if ( deltaRNew < deltaRMin ) deltaRMin = deltaRNew;
          }
          _h["dl_ph_pt"]->fill( photon.pT()/GeV );
          _h["dl_ph_eta"]->fill( photon.abseta() );
          _h["dl_ph_l_dR"]->fill( deltaRMin );
          _h["dl_l_dEta"]->fill( deltaEta(leptons[0], leptons[1]) );
          _h["dl_l_dPhi"]->fill( deltaPhi(leptons[0], leptons[1]) );
        }
      }
    }


    void finalize() {

      // Normalise histograms after the run
      for (auto &hist : _h) {
        normalize(hist.second, 1.0, true);
      }
    }


  private:
    int _mode;
    map<string, Histo1DPtr> _h;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1707015);
}
