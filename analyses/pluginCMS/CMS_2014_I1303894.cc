// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Differential cross-section of W bosons + jets in pp collisions at sqrt(s)=7 TeV
  /// @author Darin Baumgartel (darinb@cern.ch)
  ///
  /// Based on Rivet analysis originally created by Anil Singh (anil@cern.ch), Lovedeep Saini (lovedeep@cern.ch)
  class CMS_2014_I1303894 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I1303894()
      : Analysis("CMS_2014_I1303894")
    {   }


    // Book histograms and initialise projections before the run
    void init() {
      // Prompt leptons only, no test on nu flavour.
      // Projections
      const FinalState fs;
      declare(fs, "FS");

      MissingMomentum missing(fs);
      declare(missing, "MET");

      PromptFinalState pfs(fs);
      IdentifiedFinalState bareMuons(pfs);
      bareMuons.acceptIdPair(PID::MUON);
      DressedLeptons muonClusters(fs, bareMuons, -1); //, Cuts::open(), false, false);
      declare(muonClusters, "muonClusters");

      IdentifiedFinalState neutrinos(pfs);
      neutrinos.acceptIdPair(PID::NU_MU);
      declare(neutrinos, "neutrinos");

      VetoedFinalState jetFS(fs);
      jetFS.addVetoOnThisFinalState(muonClusters);
      jetFS.addVetoOnThisFinalState(neutrinos);
      jetFS.vetoNeutrinos();
      FastJets jetprojection(jetFS, FastJets::ANTIKT, 0.5);
      declare(jetprojection, "Jets");

      // Histograms
      book(_histDPhiMuJet1 ,1,1,1);
      book(_histDPhiMuJet2 ,2,1,1);
      book(_histDPhiMuJet3 ,3,1,1);
      book(_histDPhiMuJet4 ,4,1,1);

      book(_histEtaJet1 ,5,1,1);
      book(_histEtaJet2 ,6,1,1);
      book(_histEtaJet3 ,7,1,1);
      book(_histEtaJet4 ,8,1,1);

      book(_histHT1JetInc ,9,1,1);
      book(_histHT2JetInc ,10,1,1);
      book(_histHT3JetInc ,11,1,1);
      book(_histHT4JetInc ,12,1,1);

      book(_histJet30MultExc  ,13,1,1);
      book(_histJet30MultInc  ,14,1,1);

      book(_histPtJet1 ,15,1,1);
      book(_histPtJet2 ,16,1,1);
      book(_histPtJet3 ,17,1,1);
      book(_histPtJet4 ,18,1,1);

      // Counters
      book(_n_1jet, "n_1jet");
      book(_n_2jet, "n_2jet");
      book(_n_3jet, "n_3jet");
      book(_n_4jet, "n_4jet");
      book(_n_inclusivebinsummation, "n_inclusivebinsummation");
    }


    void analyze(const Event& event) {
      // Get the dressed muon
      const DressedLeptons& muonClusters = apply<DressedLeptons>(event, "muonClusters");
      int nmu = muonClusters.dressedLeptons().size();
      if (nmu < 1) vetoEvent;
      DressedLepton dressedmuon = muonClusters.dressedLeptons()[0];
      if (dressedmuon.momentum().abseta() > 2.1) vetoEvent;
      if (dressedmuon.momentum().pT() < 25.0*GeV) vetoEvent;

      // Get the muon neutrino
      //const Particles& neutrinos = apply<FinalState>(event, "neutrinos").particlesByPt();

      // Check that the muon and neutrino are not decay products of tau
      if (dressedmuon.constituentLepton().hasAncestor( PID::TAU)) vetoEvent;
      if (dressedmuon.constituentLepton().hasAncestor(-PID::TAU)) vetoEvent;

      // Get the missing momentum
      const MissingMomentum& met = apply<MissingMomentum>(event, "MET");
      const double ptmet = met.visibleMomentum().pT();
      const double phimet = (-met.visibleMomentum()).phi();

      // Calculate MET and MT(mu,MET), and remove events with MT < 50 GeV
      const double ptmuon = dressedmuon.pT();
      const double phimuon = dressedmuon.phi();
      const double mt_mumet = sqrt(2*ptmuon*ptmet*(1.0 - cos(phimet-phimuon)));

      // Remove events in MT < 50 region
      if (mt_mumet < 50*GeV) vetoEvent;

      // Loop over jets and fill pt/eta/phi quantities in vectors
      const Jets& jets_filtered = apply<FastJets>(event, "Jets").jetsByPt(0.0*GeV);
      vector<float> finaljet_pT_list, finaljet_eta_list, finaljet_phi_list;
      double htjets = 0.0;
      for (size_t ii = 0; ii < jets_filtered.size(); ++ii) {
        // Jet pT/eta/phi
        double jet_pt = jets_filtered[ii].pT();
        double jet_eta = jets_filtered[ii].eta();
        double jet_phi = jets_filtered[ii].phi();

        // Kinemetic cuts for jet acceptance
        if (fabs(jet_eta) > 2.4) continue;
        if (jet_pt < 30.0*GeV) continue;
        if (deltaR(dressedmuon, jets_filtered[ii]) < 0.5) continue;

        // Add jet to jet list and increases the HT variable
        finaljet_pT_list.push_back(jet_pt);
        finaljet_eta_list.push_back(jet_eta);
        finaljet_phi_list.push_back(jet_phi);
        htjets += fabs(jet_pt);
      }


      // Filling of histograms:
      // Fill as many jets as there are into the exclusive jet multiplicity
      if (!finaljet_pT_list.empty())
        _histJet30MultExc->fill(finaljet_pT_list.size());

      for (size_t ij = 0; ij < finaljet_pT_list.size(); ++ij) {
        _histJet30MultInc->fill(ij+1);
        _n_inclusivebinsummation->fill();
      }

      if (finaljet_pT_list.size() >= 1) {
        _histPtJet1->fill(finaljet_pT_list[0]);
        _histEtaJet1->fill(fabs(finaljet_eta_list[0]));
        _histDPhiMuJet1->fill(deltaPhi(finaljet_phi_list[0], phimuon));
        _histHT1JetInc->fill(htjets);
        _n_1jet->fill();
      }

      if (finaljet_pT_list.size() >= 2) {
        _histPtJet2->fill(finaljet_pT_list[1]);
        _histEtaJet2->fill(fabs(finaljet_eta_list[1]));
        _histDPhiMuJet2->fill(deltaPhi(finaljet_phi_list[1], phimuon));
        _histHT2JetInc->fill(htjets);
        _n_2jet->fill();
      }

      if (finaljet_pT_list.size() >= 3) {
        _histPtJet3->fill(finaljet_pT_list[2]);
        _histEtaJet3->fill(fabs(finaljet_eta_list[2]));
        _histDPhiMuJet3->fill(deltaPhi(finaljet_phi_list[2], phimuon));
        _histHT3JetInc->fill(htjets);
        _n_3jet->fill();
      }

      if (finaljet_pT_list.size() >=4 ) {
        _histPtJet4->fill(finaljet_pT_list[3]);
        _histEtaJet4->fill(fabs(finaljet_eta_list[3]));
        _histDPhiMuJet4->fill(deltaPhi(finaljet_phi_list[3], phimuon));
        _histHT4JetInc-> fill(htjets);
        _n_4jet->fill();
      }

    }


    // Finalize the histograms.
    void finalize() {

      const double inclusive_cross_section = crossSection();
      const double norm_1jet_histo = inclusive_cross_section*dbl(*_n_1jet)/sumOfWeights();
      const double norm_2jet_histo = inclusive_cross_section*dbl(*_n_2jet)/sumOfWeights();
      const double norm_3jet_histo = inclusive_cross_section*dbl(*_n_3jet)/sumOfWeights();
      const double norm_4jet_histo = inclusive_cross_section*dbl(*_n_4jet)/sumOfWeights();
      const double norm_incmultiplicity = inclusive_cross_section*dbl(*_n_inclusivebinsummation)/sumOfWeights();

      normalize(_histJet30MultExc, norm_1jet_histo);
      normalize(_histJet30MultInc, norm_incmultiplicity);

      normalize(_histPtJet1, norm_1jet_histo);
      normalize(_histHT1JetInc, norm_1jet_histo);
      normalize(_histEtaJet1, norm_1jet_histo);
      normalize(_histDPhiMuJet1, norm_1jet_histo);

      normalize(_histPtJet2, norm_2jet_histo);
      normalize(_histHT2JetInc, norm_2jet_histo);
      normalize(_histEtaJet2, norm_2jet_histo);
      normalize(_histDPhiMuJet2, norm_2jet_histo);

      normalize(_histPtJet3, norm_3jet_histo);
      normalize(_histHT3JetInc, norm_3jet_histo);
      normalize(_histEtaJet3, norm_3jet_histo);
      normalize(_histDPhiMuJet3, norm_3jet_histo);

      normalize(_histPtJet4, norm_4jet_histo);
      normalize(_histHT4JetInc, norm_4jet_histo);
      normalize(_histEtaJet4, norm_4jet_histo);
      normalize(_histDPhiMuJet4, norm_4jet_histo);
    }


  private:

    Histo1DPtr  _histJet30MultExc, _histJet30MultInc;
    Histo1DPtr  _histPtJet1, _histPtJet2, _histPtJet3, _histPtJet4;
    Histo1DPtr  _histEtaJet1, _histEtaJet2, _histEtaJet3, _histEtaJet4;
    Histo1DPtr  _histDPhiMuJet1, _histDPhiMuJet2, _histDPhiMuJet3, _histDPhiMuJet4;
    Histo1DPtr  _histHT1JetInc, _histHT2JetInc, _histHT3JetInc, _histHT4JetInc;

    CounterPtr _n_1jet, _n_2jet, _n_3jet, _n_4jet, _n_inclusivebinsummation;

  };


  DECLARE_RIVET_PLUGIN(CMS_2014_I1303894);

}
