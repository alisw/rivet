// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/XConePlugin.hh"

namespace Rivet {


  /// @brief Measurement of the jet mass for boosted top quarks at 13 TeV
  class CMS_2019_I1764472 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2019_I1764472);


    /// @name Analysis methods
    //@{

    void init() {

      // Prompt leptons
      ChargedLeptons charged_leptons;
      PromptFinalState prompt_leptons(charged_leptons);
      declare(prompt_leptons, "PromptLeptons");

      // Final state particles for jet clustering
      VetoedFinalState fs_jets;
      fs_jets.vetoNeutrinos();

      // First XCone jet clustering step
      fastjet::contrib::PseudoXConePlugin* plugin_xcone = new fastjet::contrib::PseudoXConePlugin(2, 1.2, 2.0);
      declare(FastJets(fs_jets, plugin_xcone), "FatJets");

      // Partonic tops for decay channel definition
      declare(PartonicTops(PartonicTops::DecayMode::E_MU, false), "LeptonicTops");
      declare(PartonicTops(PartonicTops::DecayMode::HADRONIC), "HadronicTops");

      // Book histograms
      book(_hist_mass, "d01-x01-y01");
      book(_hist_mass_norm, "d02-x01-y01");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Decay mode check
      const Particles& leptonicTops = apply<PartonicTops>(event, "LeptonicTops").particlesByPt();
      const Particles& hadronicTops = apply<PartonicTops>(event, "HadronicTops").particlesByPt();
      if (leptonicTops.size() != 1 || hadronicTops.size() != 1) vetoEvent;


      // Get prompt leptons
      const PromptFinalState& prompt_leptons = apply<PromptFinalState>(event, "PromptLeptons");
      const Particles & leptons = prompt_leptons.particles();
      if(leptons.empty()) vetoEvent;

      // Select leading lepton
      Particle lepton;
      for(const Particle& l : leptons){
        if(l.pT() > lepton.pT()) lepton = l;
      }
      if(lepton.pT() < 60*GeV) vetoEvent;

      // Get the fat jets
      const Jets& fatjets = applyProjection<FastJets>(event, "FatJets").jets();

      // Get index of hadronic jet by distance to lepton
      int ihad = 0;
      int ilep = 1;

      double dR0 = deltaR(lepton, fatjets.at(0));
      double dR1 = deltaR(lepton, fatjets.at(1));

      if(dR0 < dR1){
        ihad = 1;
        ilep = 0;
      }

      // Get jet constituents
      const Particles & phad = fatjets.at(ihad).particles();
      const Particles & plep = fatjets.at(ilep).particles();

      // Cluster subjets
      FinalState fs_dummy;
      fastjet::JetDefinition::Plugin* plugin_subhad = new fastjet::contrib::PseudoXConePlugin(3, 0.4, 2.0);
      fastjet::contrib::PseudoXConePlugin* plugin_sublep = new fastjet::contrib::PseudoXConePlugin(3, 0.4, 2.0);
      FastJets hadsubcluster(fs_dummy, plugin_subhad);
      FastJets lepsubcluster(fs_dummy, plugin_sublep);
      hadsubcluster.calc(phad);
      lepsubcluster.calc(plep);

      Jets subjets_had = hadsubcluster.jets();
      Jets subjets_lep = lepsubcluster.jets();

      // Subtract the lepton four vector from closest subjet if dR<0.4
      Jets subjets_had_clean;
      double dRmin_had = 0.4;
      unsigned int i_dRmin_had = 0;
      bool found_match_had = false;
      for(unsigned int i=0; i<subjets_had.size(); i++){
        double dR = deltaR(subjets_had[i], lepton);
        if(dR < dRmin_had){
          dRmin_had = dR;
          i_dRmin_had = i;
          found_match_had = true;
        }
      }
      for(unsigned int i=0; i<subjets_had.size(); i++){
        Jet subjet = subjets_had[i];
        if(found_match_had && i == i_dRmin_had) subjet = Jet(subjets_had[i].momentum()-lepton.momentum(), subjets_had[i].particles(), subjets_had[i].tags());
        subjets_had_clean.push_back(subjet);
      }
      std::sort(subjets_had_clean.begin(), subjets_had_clean.end(), cmpMomByPt);

      // do the same for lep jets
      Jets subjets_lep_clean;
      double dRmin_lep = 0.4;
      unsigned int i_dRmin_lep = 0;
      bool found_match_lep = false;
      for(unsigned int i=0; i<subjets_lep.size(); i++){
        double dR = deltaR(subjets_lep[i], lepton);
        if(dR < dRmin_lep){
          dRmin_lep = dR;
          i_dRmin_lep = i;
          found_match_lep = true;
        }
      }
      for(unsigned int i=0; i<subjets_lep.size(); i++){
        Jet subjet = subjets_lep[i];
        if(found_match_lep && i == i_dRmin_lep) subjet = Jet(subjets_lep[i].momentum()-lepton.momentum(), subjets_lep[i].particles(), subjets_lep[i].tags());
        subjets_lep_clean.push_back(subjet);
      }
      std::sort(subjets_lep_clean.begin(), subjets_lep_clean.end(), cmpMomByPt);

      // Subjet cuts
      if(subjets_had_clean.size() != 3) vetoEvent;
      if(subjets_lep_clean.size() != 3) vetoEvent;
      for (Jet jet : subjets_had_clean){
          if(jet.pT() < 30*GeV) vetoEvent;
          if(jet.abseta() > 2.5) vetoEvent;
      }

      // Combine subjets to final jets
      FourMomentum hadjet;
      for(Jet subjet : subjets_had_clean){
        if(subjet.abseta() < 2.5) hadjet += subjet.momentum();
      }
      FourMomentum lepjet;
      for(Jet subjet : subjets_lep_clean){
        if(subjet.abseta() < 2.5) lepjet += subjet.momentum();
      }

      // Jet pT cuts
      if(hadjet.pT() < 400*GeV) vetoEvent;
      if(lepjet.pT() < 10*GeV) vetoEvent;

      // m(hadjet) > m(lepjet+lepton)
      FourMomentum secondJetLepton = lepjet + lepton.momentum();
      if(hadjet.mass() < secondJetLepton.mass()) vetoEvent;

      // Fill histograms
      _hist_mass->fill(hadjet.mass()/GeV);
      _hist_mass_norm->fill(hadjet.mass()/GeV);

    }

    /// Normalise and scale histograms
    void finalize() {
      const double sf = crossSection() / femtobarn / sumOfWeights();
      scale(_hist_mass, sf);
      normalize(_hist_mass_norm, 1.0, false);
    }

    //@}


  private:

    // Histograms
    Histo1DPtr _hist_mass, _hist_mass_norm;

  };



  DECLARE_RIVET_PLUGIN(CMS_2019_I1764472);

}
