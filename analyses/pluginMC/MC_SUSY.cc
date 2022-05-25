// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {

  


  /// @brief MC validation analysis for SUSY events
  /// @author Andy Buckley
  class MC_SUSY : public Analysis {
  public:

    /// Constructor
    MC_SUSY()
      : Analysis("MC_SUSY")
    {    }


    /// @name Analysis methods
    //@{

    // Book histograms
    void init() {
      // Basic final state
      const FinalState fs((Cuts::etaIn(-4.0, 4.0) && Cuts::pT >=  10*GeV));

      // Tracks and jets
      declare(ChargedFinalState(fs), "Tracks");
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "Jets");

      IdentifiedFinalState photonfs(fs);
      photonfs.acceptId(PID::PHOTON);
      declare(photonfs, "AllPhotons");

      IdentifiedFinalState efs(fs);
      efs.acceptIdPair(PID::ELECTRON);
      declare(efs, "Electrons");

      IdentifiedFinalState mufs(fs);
      mufs.acceptIdPair(PID::MUON);
      declare(mufs, "Muons");

      MissingMomentum missing(fs);
      declare(missing, "MET");

      LeadingParticlesFinalState lpfs(fs);
      lpfs.addParticleIdPair(PID::ELECTRON);
      lpfs.addParticleIdPair(PID::MUON);
      declare(lpfs, "LeadingParticles");

      book(_hist_n_trk   ,"n-trk", 50, 0.5, 300.5);
      book(_hist_phi_trk ,"phi-trk", 50, -PI, PI);
      book(_hist_eta_trk ,"eta-trk", 50, -4, 4);
      book(_hist_pt_trk  ,"pt-trk", 100, 0.0, 1500);

      book(_hist_n_jet   ,"n-jet", 21, -0.5, 20.5);
      book(_hist_phi_jet ,"phi-jet", 50, -PI, PI);
      book(_hist_eta_jet ,"eta-jet", 50, -4, 4);
      book(_hist_pt_jet  ,"pt-jet", 100, 0.0, 1500);

      book(_hist_n_e   ,"n-e", 11, -0.5, 10.5);
      book(_hist_phi_e ,"phi-e", 50, -PI, PI);
      book(_hist_eta_e ,"eta-e", 50, -4, 4);
      book(_hist_pt_e  ,"pt-e", 100, 0.0, 500);

      book(_hist_n_mu   ,"n-mu", 11, -0.5, 10.5);
      book(_hist_phi_mu ,"phi-mu", 50, -PI, PI);
      book(_hist_eta_mu ,"eta-mu", 50, -4, 4);
      book(_hist_pt_mu  ,"pt-mu", 100, 0.0, 500);

      book(_hist_n_gamma   ,"n-gamma", 11, -0.5, 10.5);
      book(_hist_phi_gamma ,"phi-gamma", 50, -PI, PI);
      book(_hist_eta_gamma ,"eta-gamma", 50, -4, 4);
      book(_hist_pt_gamma  ,"pt-gamma", 100, 0.0, 500);

      book(_hist_n_gammaiso   ,"n-gamma-iso", 11, -0.5, 10.5);
      book(_hist_phi_gammaiso ,"phi-gamma-iso", 50, -PI, PI);
      book(_hist_eta_gammaiso ,"eta-gamma-iso", 50, -4, 4);
      book(_hist_pt_gammaiso  ,"pt-gamma-iso", 100, 0.0, 500);

      book(_hist_met ,"Etmiss", 100, 0.0, 1500);

      book(_hist_mll_ossf_ee   ,"mll-ossf-ee", 50, 0.0, 500);
      book(_hist_mll_ossf_mumu ,"mll-ossf-mumu", 50, 0.0, 500);
      book(_hist_mll_osof_emu  ,"mll-osof-emu", 50, 0.0, 500);

      book(_hist_mll_all_ossf_ee   ,"mll-all-ossf-ee", 50, 0.0, 500);
      book(_hist_mll_all_ossf_mumu ,"mll-all-ossf-mumu", 50, 0.0, 500);
      book(_hist_mll_all_osof_emu  ,"mll-all-osof-emu", 50, 0.0, 500);

      book(_hist_mll_2_ossf_ee   ,"mll-2-ossf-ee", 50, 0.0, 500);
      book(_hist_mll_2_ossf_mumu ,"mll-2-ossf-mumu", 50, 0.0, 500);
      book(_hist_mll_2_osof_emu  ,"mll-2-osof-emu", 50, 0.0, 500);

      /// @todo LSP eta, pT, phi, mass: no reliable cross-scenario LSP PID but
      /// maybe plot for all of chi^0_1, gravitino, sneutrino, gluino, ... or
      /// identify the LSP as any PID::isSUSY (?) particle with status = 1?
    }


    // Do the analysis
    void analyze(const Event& evt) {
      const FinalState& tracks = apply<FinalState>(evt, "Tracks");
      if (tracks.particles().empty()) {
        MSG_DEBUG("Failed multiplicity cut");
        vetoEvent;
      }

      // Fill track histos
      _hist_n_trk->fill(tracks.size());
      for (const Particle& t : tracks.particles()) {
        const FourMomentum& p = t.momentum();
        _hist_phi_trk->fill(mapAngleMPiToPi(p.phi()));
        _hist_eta_trk->fill(p.eta());
        _hist_pt_trk->fill(p.pT()/GeV);
      }

      // Get jets and fill jet histos
      const FastJets& jetpro = apply<FastJets>(evt, "Jets");
      const Jets jets = jetpro.jetsByPt();
      MSG_DEBUG("Jet multiplicity = " << jets.size());
      _hist_n_jet->fill(jets.size());
      for (const Jet& j : jets) {
        const FourMomentum& pj = j.momentum();
        _hist_phi_jet->fill(mapAngleMPiToPi(pj.phi()));
        _hist_eta_jet->fill(pj.eta());
        _hist_pt_jet->fill(pj.pT()/GeV);
      }

      /// @todo Resum photons around electrons

      // Fill final state electron/positron histos
      const FinalState& efs = apply<FinalState>(evt, "Electrons");
      _hist_n_e->fill(efs.size());
      vector<FourMomentum> epluses, eminuses;
      for (const Particle& e : efs.particles()) {
        const FourMomentum& p = e.momentum();
        _hist_phi_e->fill(mapAngleMPiToPi(p.phi()));
        _hist_eta_e->fill(p.eta());
        _hist_pt_e->fill(p.pT()/GeV);
        // Add sufficiently hard leptons to collections for m_ll histo
        if (p.pT()/GeV > 20) {
          if (PID::charge3(e.pid()) > 0) epluses += p; else eminuses += p;
        }
      }

      /// @todo Resum photons around muons

      // Fill final state muon/antimuon histos
      const FinalState& mufs = apply<FinalState>(evt, "Muons");
      _hist_n_mu->fill(mufs.size());
      vector<FourMomentum> mupluses, muminuses;
      for (const Particle& mu : mufs.particles()) {
        const FourMomentum& p = mu.momentum();
        _hist_phi_mu->fill(mapAngleMPiToPi(p.phi()));
        _hist_eta_mu->fill(p.eta());
        _hist_pt_mu->fill(p.pT()/GeV);
        // Add sufficiently hard leptons to collections for m_ll histo
        if (p.pT()/GeV > 20) {
          if (PID::charge3(mu.pid()) > 0) mupluses += p; else muminuses += p;
        }
      }

      // Fill final state non-isolated photon histos
      const FinalState& allphotonfs = apply<FinalState>(evt, "AllPhotons");
      _hist_n_gamma->fill(allphotonfs.size());
      Particles isolatedphotons;
      for (const Particle& ph : allphotonfs.particles()) {
        const FourMomentum& p = ph.momentum();
        _hist_phi_gamma->fill(mapAngleMPiToPi(p.phi()));
        _hist_eta_gamma->fill(p.eta());
        _hist_pt_gamma->fill(p.pT()/GeV);
        // Select isolated photons
        bool isolated = true;
        for (const Jet& j : jets) {
          if (deltaR(j.momentum(), p) < 0.2) {
            isolated = false;
            break;
          }
        }
        if (isolated) isolatedphotons += ph;
      }


      // Fill final state isolated photon histos
      _hist_n_gammaiso->fill(isolatedphotons.size());
      for (const Particle& ph_iso : isolatedphotons) {
        const FourMomentum& p = ph_iso.momentum();
        _hist_phi_gammaiso->fill(mapAngleMPiToPi(p.phi()));
        _hist_eta_gammaiso->fill(p.eta());
        _hist_pt_gammaiso->fill(p.pT()/GeV);
      }

      // Calculate and fill missing Et histos
      const MissingMomentum& met = apply<MissingMomentum>(evt, "MET");
      _hist_met->fill(met.vectorEt().mod()/GeV);

      // Choose highest-pT leptons of each sign and flavour for dilepton mass edges
      const FinalState& lpfs = apply<FinalState>(evt, "LeadingParticles");
      bool eplus_ok(false), eminus_ok(false), muplus_ok(false), muminus_ok(false);
      FourMomentum peplus, peminus, pmuplus, pmuminus;
      for (const Particle& p : lpfs.particles()) {
        // Only use leptons above 20 GeV
        if (p.pT()/GeV < 20) continue;
        // Identify the PID
        const PdgId pid = p.pid();
        if (pid == PID::ELECTRON) {
          eminus_ok = true;
          peminus = p.momentum();
        } else if (pid == PID::POSITRON) {
          eplus_ok = true;
          peplus = p.momentum();
        } else if (pid == PID::MUON) {
          muminus_ok = true;
          pmuminus = p.momentum();
        } else if (pid == PID::ANTIMUON) {
          muplus_ok = true;
          pmuplus = p.momentum();
        } else {
          throw Error("Unexpected particle type in leading particles FS!");
        }
      }
      // m_ee
      if (eminus_ok && eplus_ok) {
        const double m_ee = FourMomentum(peplus + peminus).mass();
        _hist_mll_ossf_ee->fill(m_ee/GeV);
        if (epluses.size() == 1 && eminuses.size() == 1)
          _hist_mll_2_ossf_ee->fill(m_ee/GeV);
      }
      // m_mumu
      if (muminus_ok && muplus_ok) {
        const double m_mumu = FourMomentum(pmuplus + pmuminus).mass();
        _hist_mll_ossf_mumu->fill(m_mumu/GeV);
        if (mupluses.size() == 1 && muminuses.size() == 1)
          _hist_mll_2_ossf_mumu->fill(m_mumu/GeV);
      }
      // m_emu (both configurations)
      if (eminus_ok && muplus_ok) {
        const double m_emu = FourMomentum(pmuplus + peminus).mass();
        _hist_mll_osof_emu->fill(m_emu/GeV);
        if (mupluses.size() == 1 && eminuses.size() == 1)
          _hist_mll_2_osof_emu->fill(m_emu/GeV);

      }
      if (muminus_ok && eplus_ok) {
        const double m_mue = FourMomentum(peplus + pmuminus).mass();
        _hist_mll_osof_emu->fill(m_mue/GeV);
        if (epluses.size() == 1 && muminuses.size() == 1)
          _hist_mll_2_osof_emu->fill(m_mue/GeV);
      }


      // m_ll plots using *all* electrons, positrons, muons and antimuons
      // m_ee
      for (const FourMomentum& peplus : epluses) {
        for (const FourMomentum& peminus : eminuses) {
          const double m_ee = FourMomentum(peplus + peminus).mass();
          _hist_mll_all_ossf_ee->fill(m_ee/GeV);
        }
      }
      // m_mumu
      for (const FourMomentum& pmuplus : mupluses) {
        for (const FourMomentum& pmuminus : muminuses) {
          const double m_mumu = FourMomentum(pmuplus + pmuminus).mass();
          _hist_mll_all_ossf_mumu->fill(m_mumu/GeV);
        }
      }
      // m_emu (both configurations)
      for (const FourMomentum& pmuplus : mupluses) {
        for (const FourMomentum& peminus : eminuses) {
          const double m_emu = FourMomentum(pmuplus + peminus).mass();
          _hist_mll_all_osof_emu->fill(m_emu/GeV);
        }
      }
      for (const FourMomentum& peplus : epluses) {
        for (const FourMomentum& pmuminus : muminuses) {
          const double m_mue = FourMomentum(peplus + pmuminus).mass();
          _hist_mll_all_osof_emu->fill(m_mue/GeV);
        }
      }

    }


    void finalize() {
      /// @todo Normalisations
    }

    //@}


  private:

    Histo1DPtr _hist_n_trk, _hist_phi_trk, _hist_eta_trk, _hist_pt_trk;
    Histo1DPtr _hist_n_jet, _hist_phi_jet, _hist_eta_jet, _hist_pt_jet;
    Histo1DPtr _hist_n_e, _hist_phi_e, _hist_eta_e, _hist_pt_e;
    Histo1DPtr _hist_n_mu, _hist_phi_mu, _hist_eta_mu, _hist_pt_mu;
    Histo1DPtr _hist_n_gamma, _hist_phi_gamma, _hist_eta_gamma, _hist_pt_gamma;
    Histo1DPtr _hist_n_gammaiso, _hist_phi_gammaiso, _hist_eta_gammaiso, _hist_pt_gammaiso;
    Histo1DPtr _hist_met;
    Histo1DPtr _hist_mll_2_ossf_ee, _hist_mll_2_ossf_mumu, _hist_mll_2_osof_emu;
    Histo1DPtr _hist_mll_ossf_ee, _hist_mll_ossf_mumu, _hist_mll_osof_emu;
    Histo1DPtr _hist_mll_all_ossf_ee, _hist_mll_all_ossf_mumu, _hist_mll_all_osof_emu;
  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_SUSY);

}
