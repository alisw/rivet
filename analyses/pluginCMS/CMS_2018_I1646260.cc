// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/IndirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/Smearing.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// CMS 2 soft lepton + MET in 36/fb of 13 TeV pp
  class CMS_2018_I1646260 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2018_I1646260);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON);
      SmearedParticles recoelectrons(electrons, [](const Particle& e) -> double {
          static const vector<double> ptedges = { 5., 10., 15., 20., 25., 30. };
          static const vector<double> etaedges = { 0.0, 0.8, 1.442, 1.556, 2.0, 2.5 };
          static const vector<vector<double>> effs = { { 0.336, 0.344, 0.233, 0.309, 0.243 },
                                                       { 0.412, 0.402, 0.229, 0.359, 0.287 },
                                                       { 0.465, 0.448, 0.250, 0.394, 0.327 },
                                                       { 0.496, 0.476, 0.261, 0.408, 0.341 },
                                                       { 0.503, 0.482, 0.255, 0.418, 0.352 } };
          const int ipt = binIndex(e.pT()/GeV, ptedges);
          const int ieta = binIndex(e.abseta(), etaedges);
          if (ipt < 0 || ieta < 0) return 0;
          return effs[ipt][ieta];
        }, ELECTRON_SMEAR_CMS_RUN2);
      declare(recoelectrons, "Electrons");

      PromptFinalState muons(Cuts::abspid == PID::MUON);
      SmearedParticles recomuons(muons, [](const Particle& m) -> double {
          static const vector<double> ptedges = { 3.5, 10., 15., 20., 25., 30. };
          static const vector<double> etaedges = { 0.0, 0.9, 1.2, 2.1, 2.4 };
          static const vector<vector<double>> effs = { { 0.647, 0.627, 0.610, 0.566 },
                                                       { 0.718, 0.662, 0.660, 0.629 },
                                                       { 0.739, 0.694, 0.678, 0.655 },
                                                       { 0.760, 0.725, 0.685, 0.670 },
                                                       { 0.763, 0.733, 0.723, 0.696 } };
          const int ipt = binIndex(m.pT()/GeV, ptedges);
          const int ieta = binIndex(m.abseta(), etaedges);
          if (ipt < 0 || ieta < 0) return 0;
          return effs[ipt][ieta];
        }, MUON_SMEAR_CMS_RUN2);
      declare(recomuons, "Muons");

      TauFinder taus(TauFinder::DecayMode::LEPTONIC);
      declare(taus, "Taus");

      FastJets jets4(IndirectFinalState(Cuts::abseta < 4.9), FastJets::ANTIKT, 0.4);
      SmearedJets recojets4(jets4, JET_SMEAR_CMS_RUN2, JET_BTAG_EFFS(0.8, 0.1, 0.4));
      declare(recojets4, "Jets");

      // MissingMomentum met(FinalState(Cuts::abseta < 4.9));
      // SmearedMET recomet(met, MET_SMEAR_CMS_RUN2);
      // declare(recomet, "MET");


      // Book SR counters
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 4; ++j)
          book(_srcounts_ewino[i][j], "sr_ewino_" + toString(i) + "_" + toString(j));
        for (size_t j = 0; j < 3; ++j)
          book(_srcounts_stop[i][j],   "sr_stop_" + toString(i) + "_" + toString(j));
      }

      // Cut-flow setup
      const strings cfnames = {"2mu", "mu+mu-", "pTmumu > 3 GeV", "Mmumu in [4,50] GeV",
                               "Mmumu not in [9,10.5] GeV", "pTmiss in [125, 200] GeV",
                               "mu+pTmiss trigger", "ISR jet", "HT > 100 GeV",
                               "pTmiss/HT in [0.6, 1.4]", "b-tag veto", "Mtautau veto"};
      _cutflows.addCutflow("EW", cfnames+strings{"MT < 70 GeV"});
      _cutflows.addCutflow("St", cfnames);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _cutflows.fillinit();

      // Leptons
      const Particles elecs = apply<ParticleFinder>(event, "Electrons").particlesByPt(Cuts::pT > 5*GeV && Cuts::abseta < 2.5);
      const Particles muons = apply<ParticleFinder>(event, "Muons").particlesByPt(Cuts::pT > 5*GeV && Cuts::pT < 30*GeV && Cuts::abseta < 2.4);
      const Particles leptons = sortByPt(elecs + muons);
      if (leptons.empty()) vetoEvent;
      _nevtMu += 1;

      if (leptons.size() != 2) vetoEvent;
      _cutflows.fill(1);
      if (leptons[0].charge() * leptons[1].charge() >= 0) vetoEvent;
      _cutflows.fill(2);

      // Dilepton cuts
      const FourMomentum pll = leptons[0].mom() + leptons[1].mom();
      if (pll.pT() < 3*GeV) vetoEvent;
      _cutflows.fill(3);
      const bool sameflav = (leptons[0].abspid() == leptons[1].abspid());
      if (sameflav) {
        if (!inRange(pll.mass()/GeV, 4, 50)) vetoEvent;
        _cutflows.fill(4);
        if (inRange(pll.mass()/GeV, 9, 10.5)) vetoEvent;
        _cutflows.fill(5);
      } else {
        _cutflows.fill(4);
        _cutflows.fill(5);
      }

      // Jets
      Jets jets = apply<SmearedJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.4);

      // MET
      FourMomentum p4miss;
      for (const Jet& j : jets) p4miss -= j;
      for (const Particle& e : elecs) p4miss -= e;
      const double ptmiss = p4miss.pT();
      for (const Particle& m : muons) p4miss -= m;
      const double ptmiss_mu = p4miss.pT();
      const bool mumu = (muons.size() == 2);
      if (ptmiss < (mumu ? 125*GeV : 200*GeV)) vetoEvent;
      if (ptmiss_mu < (mumu ? 125*GeV : 200*GeV)) vetoEvent;
      _cutflows.fill(6);

      // mu+pTmiss trigger (65% efficient in low-ETmiss region)
      double triggerSF = 1.0;
      if (mumu && ptmiss < 200*GeV) triggerSF = 0.65;
      _cutflows.fill(7, triggerSF);

      // ISR jet
      if (jets.empty()) vetoEvent;
      _cutflows.fill(8, triggerSF);

      // MET/HT
      const double ht = sum(jets, Kin::pT, 0.0);
      if (ht < 100*GeV) vetoEvent;
      _cutflows.fill(9, triggerSF);
      if (!inRange(ptmiss/ht, 0.6, 1.4)) vetoEvent;
      _cutflows.fill(10, triggerSF);

      // b-jet veto
      if (any(jets, hasBTag(Cuts::pT > 5*GeV))) vetoEvent; //< b-jet veto with ad hoc tagging threshold
      _cutflows.fill(11, triggerSF);

      // Tau veto
      const Particles taus = apply<ParticleFinder>(event, "Taus").particlesByPt();
      if (taus.size() >= 2) {
        const double mtt = (taus[0].mom() + taus[1].mom()).mass();
        if (mtt < 160*GeV) vetoEvent;
      }
      _cutflows.fill(12, triggerSF);

      // EWino SR (ee, mumu)
      if (sameflav) {
        for (const Particle& l : leptons) {
          if (l.abspid() == PID::MUON and l.pT() < 5*GeV) vetoEvent;
          if (mT(l.mom(), p4miss) > 70*GeV) vetoEvent;
        }
        _cutflows["EW"].fill(13, triggerSF);
        static const vector<double> ptmissedges_ewino = {125., 200., 250., DBL_MAX};
        static const vector<double> mlledges_ewino = {4., 9., 10.5, 20., 30., 50.};
        const int iptm = binIndex(ptmiss/GeV, ptmissedges_ewino);
        const int imll = binIndex(pll.mass()/GeV, mlledges_ewino);
        _srcounts_ewino[iptm][imll < 1 ? 0 : imll-1]->fill(triggerSF);

      } else {

        // Stop SR
        if (leptons[0].abspid() == PID::MUON and leptons[0].pT() < 5*GeV) vetoEvent;
        _cutflows["St"].fill(13, triggerSF);
        static const vector<double> ptmissedges_stop = {125., 200., 300., DBL_MAX};
        static const vector<double> ptledges_stop = {5., 12., 20., 30.};
        const int iptm = binIndex(ptmiss/GeV, ptmissedges_stop);
        const int iptl = binIndex(leptons[0].pT()/GeV, ptledges_stop);
        _srcounts_stop[iptm][iptl]->fill(triggerSF);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      MSG_INFO("Num events with >= 1 muon = " << _nevtMu << " / " << numEvents());

      const double sf = 35.9*crossSection()/femtobarn/sumOfWeights();
      for (size_t i = 0; i < 3; ++i) {
        scale(_srcounts_ewino[i], sf);
        scale(_srcounts_stop[i], sf);
      }
      _cutflows.scale(sf);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflows);
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _srcounts_ewino[3][4], _srcounts_stop[3][3];
    //@}

    /// Cut-flows
    int _nevtMu = 0;
    Cutflows _cutflows;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_I1646260);


}
