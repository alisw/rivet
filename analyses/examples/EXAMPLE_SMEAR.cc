// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedMET.hh"

namespace Rivet {


  class EXAMPLE_SMEAR : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_SMEAR);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      MissingMomentum mm(Cuts::abseta < 5);
      declare(mm, "MET0");

      SmearedMET smm1(mm, MET_SMEAR_ATLAS_RUN2);
      declare(smm1, "MET1");

      SmearedMET smm2(mm, [](const Vector3& met, double){ return P3_SMEAR_LEN_GAUSS(met, 0.1*met.mod()); });
      declare(smm2, "MET2");


      FastJets fj(FinalState(Cuts::abseta < 5), FastJets::ANTIKT, 0.4);
      declare(fj, "Jets0");

      SmearedJets sj1(fj, JET_SMEAR_IDENTITY);
      declare(sj1, "Jets1");

      SmearedJets sj2(fj, JET_SMEAR_ATLAS_RUN2,
                      [](const Jet& j){ return j.bTagged() ? 0.7*(1 - exp(-j.pT()/(10*GeV))) : 0.01; } );
      declare(sj2, "Jets2");

      SmearedJets sj3(fj,
                      JET_SMEAR_CMS_RUN2,
                      JET_BTAG_EFFS(0.7, 0.1, 0.01),
                      JET_CTAG_PERFECT,
                      JET_EFF_CONST(0.8));
      declare(sj3, "Jets3");


      IdentifiedFinalState photons(Cuts::abseta < 5, PID::PHOTON);


      IdentifiedFinalState truthelectrons(Cuts::abseta < 5 && Cuts::pT > 10*GeV, {{PID::ELECTRON, PID::POSITRON}});
      declare(truthelectrons, "Electrons0");
      DressedLeptons dressedelectrons(photons, truthelectrons, 0.2);
      declare(dressedelectrons, "Electrons1");
      SmearedParticles recoelectrons(dressedelectrons, ELECTRON_RECOEFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2);
      declare(recoelectrons, "Electrons2");

      IdentifiedFinalState truthmuons(Cuts::abseta < 5 && Cuts::pT > 10*GeV, {{PID::MUON, PID::ANTIMUON}});
      declare(truthmuons, "Muons0");
      DressedLeptons dressedmuons(photons, truthmuons, 0.2);
      declare(dressedmuons, "Muons1");
      SmearedParticles recomuons(dressedmuons, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2);
      declare(recomuons, "Muons2");

      TauFinder truthtaus(TauFinder::DecayMode::ANY, Cuts::abseta < 5 && Cuts::pT > 10*GeV);
      declare(truthtaus, "Taus0");
      DressedLeptons dressedtaus(photons, truthtaus, 0.2);
      declare(dressedtaus, "Taus1");
      SmearedParticles recotaus(dressedtaus, TAU_EFF_ATLAS_RUN2, TAU_SMEAR_ATLAS_RUN2);
      declare(recotaus, "Taus2");


      book(_h_met_true ,"met_true", 30, 0.0, 120);
      book(_h_met_reco ,"met_reco", 30, 0.0, 120);

      book(_h_nj_true ,"jet_N_true", 10, -0.5, 9.5);
      book(_h_nj_reco ,"jet_N_reco", 10, -0.5, 9.5);
      book(_h_j1pt_true ,"jet_pt1_true", 30, 0.0, 120);
      book(_h_j1pt_reco ,"jet_pt1_reco", 30, 0.0, 120);
      book(_h_j1eta_true ,"jet_eta1_true", 20, -5.0, 5.0);
      book(_h_j1eta_reco ,"jet_eta1_reco", 20, -5.0, 5.0);

      book(_h_ne_true ,"elec_N_true", 5, -0.5, 4.5);
      book(_h_ne_reco ,"elec_N_reco", 5, -0.5, 4.5);
      book(_h_e1pt_true ,"elec_pt1_true", 30, 0, 120);
      book(_h_e1pt_reco ,"elec_pt1_reco", 30, 0, 120);
      book(_h_e1eta_true ,"elec_eta1_true", 20, -5.0, 5.0);
      book(_h_e1eta_reco ,"elec_eta1_reco", 20, -5.0, 5.0);

      book(_h_nm_true ,"muon_N_true", 5, -0.5, 4.5);
      book(_h_nm_reco ,"muon_N_reco", 5, -0.5, 4.5);
      book(_h_m1pt_true ,"muon_pt1_true", 30, 0, 120);
      book(_h_m1pt_reco ,"muon_pt1_reco", 30, 0, 120);
      book(_h_m1eta_true ,"muon_eta1_true", 20, -5.0, 5.0);
      book(_h_m1eta_reco ,"muon_eta1_reco", 20, -5.0, 5.0);

      book(_h_nt_true ,"tau_N_true", 5, -0.5, 4.5);
      book(_h_nt_reco ,"tau_N_reco", 5, -0.5, 4.5);
      book(_h_t1pt_true ,"tau_pt1_true", 30, 0, 120);
      book(_h_t1pt_reco ,"tau_pt1_reco", 30, 0, 120);
      book(_h_t1eta_true ,"tau_eta1_true", 20, -5.0, 5.0);
      book(_h_t1eta_reco ,"tau_eta1_reco", 20, -5.0, 5.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      const Vector3 met0 = apply<MissingMomentum>(event, "MET0").vectorEt();
      const Vector3 met1 = apply<SmearedMET>(event, "MET1").vectorEt();
      const Vector3 met2 = apply<SmearedMET>(event, "MET2").vectorEt();
      MSG_DEBUG("MET = " << met0.mod()/GeV << ", " << met1.mod()/GeV << ", " << met2.mod()/GeV << " GeV");
      _h_met_true->fill(met0.mod()/GeV, weight);
      _h_met_reco->fill(met1.mod()/GeV, weight);
      if (met0.perp() > 0 && met1.perp() > 0 && deltaPhi(met0, met1) > 0.1) {
        MSG_WARNING("Large MET phi change: " << met0.phi()  << " -> " << met1.phi() <<
                    "; dphi = " << deltaPhi(met0, met1));
      }

      const Jets jets0 = apply<JetAlg>(event, "Jets0").jetsByPt(Cuts::pT > 10*GeV);
      const Jets jets1 = apply<JetAlg>(event, "Jets1").jetsByPt(Cuts::pT > 10*GeV);
      const Jets jets2 = apply<JetAlg>(event, "Jets2").jetsByPt(Cuts::pT > 10*GeV);
      const Jets jets3 = apply<JetAlg>(event, "Jets3").jetsByPt(Cuts::pT > 10*GeV);
      MSG_DEBUG("Numbers of jets = " << jets0.size() << " true; "
               << jets1.size() << ", " << jets2.size() << ", " << jets3.size());
      if (!jets0.empty() && !jets2.empty() && deltaPhi(jets0[0], jets2[0]) > 0.1) {
        MSG_DEBUG("Large jet1 phi change (could be a different truth-jet): " <<
                  jets0[0].phi() << " -> " << jets2[0].phi() <<
                  "; dphi = " << deltaPhi(jets0[0], jets2[0]) <<
                  "; pT = " << jets0[0].pT()/GeV << " -> " << jets2[0].pT()/GeV);
      }
      _h_nj_true->fill(jets0.size(), weight);
      _h_nj_reco->fill(jets2.size(), weight);
      if (!jets0.empty()) {
        _h_j1pt_true->fill(jets0.front().pT()/GeV, weight);
        _h_j1eta_true->fill(jets0.front().eta(), weight);
      }
      if (!jets2.empty()) {
        _h_j1pt_reco->fill(jets2.front().pT()/GeV, weight);
        _h_j1eta_reco->fill(jets2.front().eta(), weight);
      }

      const Particles& elecs1 = apply<ParticleFinder>(event, "Electrons1").particlesByPt();
      const Particles& elecs2 = apply<ParticleFinder>(event, "Electrons2").particlesByPt();
      MSG_DEBUG("Numbers of electrons = " << elecs1.size() << " true; " << elecs2.size() << " reco");
      _h_ne_true->fill(elecs1.size(), weight);
      _h_ne_reco->fill(elecs2.size(), weight);
      if (!elecs1.empty()) {
        _h_e1pt_true->fill(elecs1.front().pT()/GeV, weight);
        _h_e1eta_true->fill(elecs1.front().eta(), weight);
      }
      if (!elecs2.empty()) {
        _h_e1pt_reco->fill(elecs2.front().pT()/GeV, weight);
        _h_e1eta_reco->fill(elecs2.front().eta(), weight);
      }

      const Particles& muons1 = apply<ParticleFinder>(event, "Muons1").particlesByPt();
      const Particles& muons2 = apply<ParticleFinder>(event, "Muons2").particlesByPt();
      MSG_DEBUG("Numbers of muons = " << muons1.size() << " true; " << muons2.size() << " reco");
      _h_nm_true->fill(muons1.size(), weight);
      _h_nm_reco->fill(muons2.size(), weight);
      if (!muons1.empty()) {
        _h_m1pt_true->fill(muons1.front().pT()/GeV, weight);
        _h_m1eta_true->fill(muons1.front().eta(), weight);
      }
      if (!muons2.empty()) {
        _h_m1pt_reco->fill(muons2.front().pT()/GeV, weight);
        _h_m1eta_reco->fill(muons2.front().eta(), weight);
      }

      const Particles& taus1 = apply<ParticleFinder>(event, "Taus1").particlesByPt();
      const Particles& taus2 = apply<ParticleFinder>(event, "Taus2").particlesByPt();
      MSG_DEBUG("Numbers of taus = " << taus1.size() << " true; " << taus2.size() << " reco");
      _h_nt_true->fill(taus1.size(), weight);
      _h_nt_reco->fill(taus2.size(), weight);
      if (!taus1.empty()) {
        _h_t1pt_true->fill(taus1.front().pT()/GeV, weight);
        _h_t1eta_true->fill(taus1.front().eta(), weight);
      }
      if (!taus2.empty()) {
        _h_t1pt_reco->fill(taus2.front().pT()/GeV, weight);
        _h_t1eta_reco->fill(taus2.front().eta(), weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_met_true);
      normalize(_h_met_reco);

      normalize(_h_nj_true);
      normalize(_h_nj_reco);
      normalize(_h_j1pt_true, 1-_h_nj_true->bin(0).area());
      normalize(_h_j1pt_reco, 1-_h_nj_reco->bin(0).area());
      normalize(_h_j1eta_true, 1-_h_nj_true->bin(0).area());
      normalize(_h_j1eta_reco, 1-_h_nj_reco->bin(0).area());

      normalize(_h_ne_true);
      normalize(_h_ne_reco);
      normalize(_h_e1pt_true, 1-_h_ne_true->bin(0).area());
      normalize(_h_e1pt_reco, 1-_h_ne_reco->bin(0).area());
      normalize(_h_e1eta_true, 1-_h_ne_true->bin(0).area());
      normalize(_h_e1eta_reco, 1-_h_ne_reco->bin(0).area());

      normalize(_h_nm_true);
      normalize(_h_nm_reco);
      normalize(_h_m1pt_true, 1-_h_nm_true->bin(0).area());
      normalize(_h_m1pt_reco, 1-_h_nm_reco->bin(0).area());
      normalize(_h_m1eta_true, 1-_h_nm_true->bin(0).area());
      normalize(_h_m1eta_reco, 1-_h_nm_reco->bin(0).area());

      normalize(_h_nt_true);
      normalize(_h_nt_reco);
      normalize(_h_t1pt_true, 1-_h_nt_true->bin(0).area());
      normalize(_h_t1pt_reco, 1-_h_nt_reco->bin(0).area());
      normalize(_h_t1eta_true, 1-_h_nt_true->bin(0).area());
      normalize(_h_t1eta_reco, 1-_h_nt_reco->bin(0).area());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_met_true, _h_met_reco;
    Histo1DPtr _h_nj_true, _h_nj_reco, _h_ne_true, _h_ne_reco,  _h_nm_true, _h_nm_reco,  _h_nt_true, _h_nt_reco;
    Histo1DPtr _h_j1pt_true, _h_j1pt_reco, _h_e1pt_true, _h_e1pt_reco,  _h_m1pt_true, _h_m1pt_reco,  _h_t1pt_true, _h_t1pt_reco;
    Histo1DPtr _h_j1eta_true, _h_j1eta_reco, _h_e1eta_true, _h_e1eta_reco,  _h_m1eta_true, _h_m1eta_reco,  _h_t1eta_true, _h_t1eta_reco;
    //@}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(EXAMPLE_SMEAR);


}
