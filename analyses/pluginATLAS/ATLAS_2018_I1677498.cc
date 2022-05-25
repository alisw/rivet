#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// WWbb at 13 TeV
  class ATLAS_2018_I1677498 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1677498);


    /// Book cuts and projections
    void init() {

      // All final state particles
      FinalState fs(Cuts::abseta < 5.0);

      PromptFinalState photons(Cuts::abspid == PID::PHOTON, true); // true accepts tau decays

      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true); // true accepts tau decays
      DressedLeptons elecs(photons, bare_el, 0.1, Cuts::pT > 7*GeV && Cuts::abseta < 2.47);
      declare(elecs, "elecs");

      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true); // accepts tau decays
      DressedLeptons muons(photons, bare_mu, 0.1, Cuts::pT > 6*GeV && Cuts::abseta < 2.5);
      declare(muons, "muons");

      FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "jets");

      book(_h, 3, 1, 1);
    }


    void analyze(const Event& event) {

      // Identify bjets (signal), light jets (OR) and soft b-jets (veto)
      size_t soft_bjets = 0;
      Jets bjets, lightjets;
      for (const Jet& jet : apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 5*GeV)) {
        bool isBjet = jet.bTagged(Cuts::pT > 5*GeV);
        if (isBjet) soft_bjets += 1;
        if (jet.abseta() < 2.5) {
          if ( isBjet && jet.pT() > 25*GeV) bjets += jet;
          if (!isBjet && jet.pT() > 20*GeV) lightjets += jet;
        }
      }
      if (soft_bjets != 2) vetoEvent;
      if (bjets.size() != 2) vetoEvent;

      // Get dressed leptons
      vector<DressedLepton> leptons;
      for (auto& lep : apply<DressedLeptons>(event, "muons").dressedLeptons()) { leptons.push_back(lep); }
      for (auto& lep : apply<DressedLeptons>(event, "elecs").dressedLeptons()) { leptons.push_back(lep); }

      // 1. Find which light jets survive OR
      for (const auto& lep : leptons) {
        ifilter_discard(lightjets, [&](const Jet& jet) {
          return deltaR(jet, lep) < 0.2 && (lep.abspid() == PID::ELECTRON || lep.pT()/jet.pT() > 0.7);
        });
      }

      // 2. Find which leptons survive the OR and apply signal selection
      for (const auto& jet : (lightjets + bjets)) {
        ifilter_discard(leptons, [&](const DressedLepton& lep) {
          return lep.pT() < 28*GeV || deltaR(jet, lep) < min(0.4, 0.04+10*GeV/lep.pT());
        });
      }

      if (leptons.size() != 2) vetoEvent;
      std::sort(leptons.begin(), leptons.end(), cmpMomByPt);

      // Z veto
      const size_t nEl = count(leptons, [](const DressedLepton& lep) { return  lep.abspid() == PID::ELECTRON; });
      const double mll = (leptons[0].mom() + leptons[1].mom()).mass();
      if (nEl != 1 && !(fabs(mll - 91*GeV) > 15*GeV && mll > 10*GeV)) vetoEvent;

      const double m00 = (leptons[0].mom() + bjets[0].mom()).mass();
      const double m10 = (leptons[1].mom() + bjets[0].mom()).mass();
      const double m01 = (leptons[0].mom() + bjets[1].mom()).mass();
      const double m11 = (leptons[1].mom() + bjets[1].mom()).mass();
      double minimax = min( max(m00,m11), max(m10,m01) );
      minimax = min(minimax, 419.); // put overflow in the last bin
      _h->fill(minimax/GeV);
    }


    /// Finalise
    void finalize() {
      normalize(_h);
    }


  private:

    /// Histogram
    Histo1DPtr _h;

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1677498);

}
