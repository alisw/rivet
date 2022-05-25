// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief b-fragmentation at 13 TeV
  class ATLAS_2021_I1913061 : public Analysis {

  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2021_I1913061);

    /// @name Analysis methods
    //@{

    void init() {

      //Jet building: anti-kT R = 0.4, including muons and neutrinos.
      FinalState photons(Cuts::abspid == PID::PHOTON);

      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true);
      DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, Cuts::abseta < 2.5, true);

      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true);
      DressedLeptons all_dressed_el(photons, bare_el, 0.1, Cuts::abseta < 2.5, true);

      VetoedFinalState vfs(FinalState(Cuts::abseta < 4.5));
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);
      
      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jets, "JETS");

      //Charged B mesons
      UnstableParticles bpm_fs(Cuts::abspid == 521);
      declare(bpm_fs, "BPM_FS");

      // Book histograms
      book(_h["zFrag_pt01"], 1,1,1);
      book(_h["ptRel_pt01"], 2,1,1);
      book(_h["zFrag_pt02"], 3,1,1);
      book(_h["ptRel_pt02"], 4,1,1);
      book(_h["zFrag_pt03"], 5,1,1);
      book(_h["ptRel_pt03"], 6,1,1);
      book(_p["zFrag"], 7,1,1);
      book(_p["ptRel"], 8,1,1);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const Jets& jets = apply<FastJets>(event, "JETS").jetsByPt(7*GeV);
      const Particles& bpmFS = apply<UnstableParticles>(event, "BPM_FS").particlesByPt();      

      //Preselect B mesons in J/psi K decay channel
      Particles goodHadrons;

      for (const Particle& p : bpmFS) {

        bool passKaon = 0;
        bool passMuon1 = 0;
        bool passMuon2 = 0;

        const Particles& bDecays = p.children();

        for (const Particle& s : bDecays){
          
          if (s.abspid() == PID::KPLUS && s.pt() > 4*GeV) passKaon = 1;
          if (s.abspid() == PID::JPSI) {

            const Particles& jpsiDecays = s.children();

            for (const Particle& m : jpsiDecays){

              if (m.pid() == PID::ANTIMUON && m.pt() > 6*GeV) passMuon1 = 1;
              if (m.pid() == PID::MUON && m.pt() > 6*GeV) passMuon2 = 1;
            }
          }
        }

        if (passKaon && passMuon1 && passMuon2) goodHadrons += p;
      }
          
      //Preselect jets passing kinematic cuts
      Jets goodJets;
      for (size_t i = 0; i < jets.size(); ++i) {
        if (jets[i].pT() <= 20*GeV) continue;
        if (jets[i].abseta() >= 2.1) continue;
        bool overlaps = false;
        for (size_t j = i + 1; j < jets.size(); ++j) {
          if (jets[j].pT() > 20*GeV && deltaR(jets[i], jets[j]) < 0.8) {
            overlaps = true; break;
          }
        }
        if (!overlaps)  goodJets += jets[i];
      }

      //Associate the jets to the B hadrons.
      for (const Jet& thisJet : goodJets) {

        Particles jetHads;
        for (const Particle& thisHad : goodHadrons) {
          if (deltaR(thisJet, thisHad) < 0.4) jetHads += thisHad;
        }

        if (jetHads.size() == 1){
          const Particle& jetHadron = jetHads[0];

          //Out-of-cone correction for B-decay products
          Vector3 jetVector = momentum3(thisJet); Vector3 jetVector0 = jetVector;
          Vector3 hadronVector = momentum3(jetHadron);
          Vector3 kaonVector; Vector3 muonVector1; Vector3 muonVector2;
          
          const Particles& bDecays = jetHadron.children();
          for (const Particle& s : bDecays){
            
            if (s.abspid() == PID::KPLUS) kaonVector = momentum3(s);
            if (s.abspid() == PID::JPSI) {

              const Particles& jpsiDecays = s.children();

              for (const Particle& m : jpsiDecays){

                if (m.pid() == PID::ANTIMUON && m.pt() > 6*GeV) muonVector1 = momentum3(m);
                if (m.pid() == PID::MUON && m.pt() > 6*GeV) muonVector2 = momentum3(m);
              }
            }
          }

          if (deltaR(muonVector1, jetVector0) > 0.4) jetVector += muonVector1;
          if (deltaR(muonVector2, jetVector0) > 0.4) jetVector += muonVector2;
          if (deltaR(kaonVector, jetVector0) > 0.4) jetVector += kaonVector;

          if (jetVector.perp() <= 30.0*GeV) vetoEvent;
          if (jetVector.abseta() >= 2.1)    vetoEvent;
          if (deltaR(jetVector, hadronVector) >= 0.4)  vetoEvent;

          //Longitudinal and transverse profiles
          double zFrag = hadronVector.dot(jetVector)/jetVector.mod2();
          double ptRel = (hadronVector.cross(jetVector)).mod()/jetVector.mod();
          
          if (jetVector.perp() >= 50*GeV && jetVector.perp() < 70*GeV) {
            if (zFrag <= 0.23) zFrag = 0.24;
            if (zFrag >= 1.00) zFrag = 0.99;
            if (ptRel >= 8.00) ptRel = 7.90; 

            _h["zFrag_pt01"]->fill(zFrag);
            _h["ptRel_pt01"]->fill(ptRel);
          }

          if (jetVector.perp() >= 70*GeV && jetVector.perp() < 100*GeV) {
            if (zFrag <= 0.23) zFrag = 0.24;
            if (zFrag >= 1.00) zFrag = 0.99;
            if (ptRel >= 10.0) ptRel = 9.90;

            _h["zFrag_pt02"]->fill(zFrag);
            _h["ptRel_pt02"]->fill(ptRel);
          }

          if (jetVector.perp() >= 100*GeV) {
            if (zFrag <= 0.16) zFrag = 0.17;
            if (zFrag >= 1.00) zFrag = 0.99;
            if (ptRel >= 14.0) ptRel = 13.9;

            _h["zFrag_pt03"]->fill(zFrag);
            _h["ptRel_pt03"]->fill(ptRel);
          }

          double ptFill = jetVector.perp()/GeV;
          if (ptFill >= 150.)  ptFill = 125.;
		
          // For pTrel, we need to weight the mean by 1/binWidth 
          // (what is being plotted is the mean of the histogram, not the variable!)
          // This doesn't apply to z since all bins have the same width.
          double binWidth = 0.0;
          double transBins[12] = { 0.0, 0.5, 1.0, 1.5, 2.2, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 14.0 };

          for (size_t iBin = 0; iBin < 11; ++iBin) {
            if (ptRel >= transBins[iBin] && ptRel < transBins[iBin+1])  binWidth = transBins[iBin+1]-transBins[iBin];
          }

          _p["zFrag"]->fill(ptFill, zFrag);
          _p["ptRel"]->fill(ptFill, ptRel, 1./binWidth);

        }
      }
    }
    
    void finalize() {
      normalize(_h);
    }
    
  private:
  
    //Histograms
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2021_I1913061);
}
