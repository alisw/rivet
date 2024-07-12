#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief same-sign WW at 13 TeV
  class ATLAS_2019_I1738841 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1738841);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState allLeps(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      PromptFinalState promptLeps(allLeps);
      FinalState photons(Cuts::abspid == PID::PHOTON);

      Cut dressedLep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 27*GeV);
      const DressedLeptons dressedLeps(photons, promptLeps, 0.1, dressedLep_cuts, true);
      declare(dressedLeps, "dressedLeptons");

      // veto on leptons from prompt tau decays
      VetoedFinalState lepsFromTaus(PromptFinalState(allLeps, true));
      lepsFromTaus.addVetoOnThisFinalState(promptLeps);
      const DressedLeptons vetoLeps(photons, lepsFromTaus, 0.1, dressedLep_cuts, true);
      declare(vetoLeps, "vetoLeptons");

      declare(MissingMomentum(), "eTmiss");

      VetoedFinalState vfs(FinalState(Cuts::abseta < 4.5));
      vfs.addVetoOnThisFinalState(dressedLeps);
      const FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::NONE);
      declare(jets, "jets");

      // histograms
      book(_c, 1, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (apply<DressedLeptons>(event, "vetoLeptons").dressedLeptons().size())  vetoEvent;

      const vector<DressedLepton> &dressedLeptons = apply<DressedLeptons>(event, "dressedLeptons").dressedLeptons();
      if (dressedLeptons.size() != 2) vetoEvent;
      const Particle& lep1 = dressedLeptons[0];
      const Particle& lep2 = dressedLeptons[1];

      const Jets& jets = apply<FastJets>(event,"jets").jetsByPt(35.);
      if (jets.size() < 2) vetoEvent;

      // pT of the two leading jets
      const Jet& tagJet1 = jets[0];
      const Jet& tagJet2 = jets[1];
      if (tagJet1.pT() < 65.)  vetoEvent;
      if (tagJet2.pT() < 35.)  vetoEvent;

      // leptons isolation from jets (dR >= 0.3)
      float f_lj_dRmin = 999;
      for (const Jet& j : jets) {
        float dR = deltaR(j, lep1);
        if(dR < f_lj_dRmin)  f_lj_dRmin = dR;
        dR = deltaR(j, lep2);
        if(dR < f_lj_dRmin)  f_lj_dRmin = dR;
      }
      if (f_lj_dRmin < 0.3) vetoEvent;

      // same sign leptons
      if (lep1.pid() * lep2.pid() < 0) vetoEvent;

      // lepton isolation
      const float f_ll_dR = deltaR(lep1.momentum(), lep2.momentum());
      if (f_ll_dR < 0.3) vetoEvent;

      // m_ll > 20 GeV
      const float f_ll_m = (lep1.momentum() + lep2.momentum()).mass()/GeV;
      if (f_ll_m < 20.) vetoEvent;

      // MET > 30 GeV
      const double ETmiss = apply<MissingMomentum>(event, "eTmiss").missingPt()/GeV;
      if (ETmiss < 30.)  vetoEvent;

      // m_jj >= 500 GeV
      const float f_tjets_m = (tagJet1.momentum() + tagJet2.momentum()).mass()/GeV;
      if (f_tjets_m < 500.)  vetoEvent;

      // dY_jj >= 2
      const float f_tjets_dY = deltaRap(tagJet1, tagJet2);
      if (f_tjets_dY < 2)  vetoEvent;

      // --- fill histograms ---
      _c->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c, crossSection() / femtobarn / sumOfWeights());
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _c;
    //@}

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1738841);
}
