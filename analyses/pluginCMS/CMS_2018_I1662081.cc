#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  class CMS_2018_I1662081 : public Analysis {
  public:

    // Minimal constructor
    CMS_2018_I1662081()
      : Analysis("CMS_2018_I1662081") {}


    // Set up projections and book histograms
    void init() {
      // Complete final state
      FinalState fs( (Cuts::abseta < 5) and (Cuts::pT > 0.0*MeV) );

      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);
      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);
      Cut looseLeptonCuts = Cuts::pt > 15*GeV && Cuts::abseta < 2.4;

      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, looseLeptonCuts, true);
      declare(dressed_leptons, "DressedLeptons");

      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      declare(FastJets(fsForJets, FastJets::ANTIKT, 0.4), "Jets");

      // Projections for MET
      declare(MissingMomentum(fs), "MET");

      // Booking of histograms
      book(_hist_norm_met , 4, 1, 1);
      book(_hist_norm_ht  , 2, 1, 1);
      book(_hist_norm_st  , 3, 1, 1);
      book(_hist_norm_wpt , 5, 1, 1);
      book(_hist_norm_njets , 1, 1, 1);
      book(_hist_norm_lpt , 6, 1, 1);
      book(_hist_norm_labseta , 7, 1, 1);

      book(_hist_abs_met , 11, 1, 1);
      book(_hist_abs_ht  , 9, 1, 1);
      book(_hist_abs_st  , 10, 1, 1);
      book(_hist_abs_wpt , 12, 1, 1);
      book(_hist_abs_njets , 8, 1, 1);
      book(_hist_abs_lpt , 13, 1, 1);
      book(_hist_abs_labseta , 14, 1, 1);

    }


    // per event analysis
    void analyze(const Event& event) {

      // Lepton veto selection
      const DressedLeptons& dressed_leptons = applyProjection<DressedLeptons>(event, "DressedLeptons");
      if (dressed_leptons.dressedLeptons().size() != 1) vetoEvent;

      // Signal lepton selection
      FourMomentum lepton = dressed_leptons.dressedLeptons()[0];

      const double leptonPt = lepton.pT();
      const double leptonAbsEta = std::abs( lepton.eta() );
      if (leptonPt <= 26*GeV || leptonAbsEta >= 2.4) vetoEvent;

      // Jet selection
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jets(Cuts::abseta < 2.4 && Cuts::pT > 20*GeV);
      Jets cleanedJets;
      unsigned int nJetsAbove30GeV = 0;
      unsigned int nJetsAbove20GeV = 0;
      unsigned int nBJetsAbove30GeV = 0;
      unsigned int nBJetsAbove20GeV = 0;
      for (const Jet& j : jets) {
        cleanedJets.push_back( j );
        ++nJetsAbove20GeV;
        if ( j.pT() > 30*GeV) ++nJetsAbove30GeV;

        if ( j.bTagged() ) {
          ++nBJetsAbove20GeV;
          if ( j.pT() > 30*GeV) ++nBJetsAbove30GeV;
        }
      }

      if ( nJetsAbove30GeV < 3 || nJetsAbove20GeV < 4 ) vetoEvent;
      if ( nBJetsAbove30GeV < 1 || nBJetsAbove20GeV < 2 ) vetoEvent;

      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _hist_norm_met->fill(met.visibleMomentum().pT()/GeV);
      _hist_abs_met->fill(met.visibleMomentum().pT()/GeV);

      // HT and ST
      double ht = 0.0;
      for (const Jet& j : cleanedJets) ht += j.pT();

      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _hist_norm_ht->fill(ht/GeV);
      _hist_norm_st->fill(st/GeV);
      _hist_abs_ht->fill(ht/GeV);
      _hist_abs_st->fill(st/GeV);

      // WPT
      FourMomentum w = lepton - met.visibleMomentum();
      _hist_norm_wpt->fill(w.pT()/GeV);
      _hist_abs_wpt->fill(w.pT()/GeV);

      // Lepton pt and eta
      _hist_norm_lpt->fill( leptonPt/GeV);
      _hist_norm_labseta->fill( leptonAbsEta/GeV);

      _hist_abs_lpt->fill( leptonPt/GeV);
      _hist_abs_labseta->fill( leptonAbsEta/GeV);

      // NJets
      _hist_norm_njets->fill( cleanedJets.size());
      _hist_abs_njets->fill( cleanedJets.size());

    }


    void finalize() {
      normalize(_hist_norm_met);
      normalize(_hist_norm_ht);
      normalize(_hist_norm_st);
      normalize(_hist_norm_wpt);
      normalize(_hist_norm_njets);
      normalize(_hist_norm_lpt);
      normalize(_hist_norm_labseta);

      scale(_hist_abs_met, crossSection()/picobarn / sumOfWeights());
      scale(_hist_abs_ht, crossSection()/picobarn / sumOfWeights());
      scale(_hist_abs_st, crossSection()/picobarn / sumOfWeights());
      scale(_hist_abs_wpt, crossSection()/picobarn / sumOfWeights());
      scale(_hist_abs_njets, crossSection()/picobarn / sumOfWeights());
      scale(_hist_abs_lpt, crossSection()/picobarn / sumOfWeights());
      scale(_hist_abs_labseta, crossSection()/picobarn / sumOfWeights());

    }


  private:
    Histo1DPtr _hist_norm_met, _hist_norm_ht, _hist_norm_st, _hist_norm_wpt, _hist_norm_njets, _hist_norm_lpt, _hist_norm_labseta;
    Histo1DPtr _hist_abs_met, _hist_abs_ht, _hist_abs_st, _hist_abs_wpt, _hist_abs_njets, _hist_abs_lpt, _hist_abs_labseta;

  };


  RIVET_DECLARE_PLUGIN(CMS_2018_I1662081);

}
