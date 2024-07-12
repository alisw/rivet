// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "TFile.h"
#include "TTree.h"

namespace Rivet {


  /// @brief An example of writing a ROOT file with ntuple data per event
  class EXAMPLE_NTUPLE_ROOT : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_NTUPLE_ROOT);


    /// @name Analysis methods
    /// @{

    /// Initialise projections and output before the run
    void init() {

      // Get jets, leptons, MET
      const FinalState fs(Cuts::abseta < 4.9);
      declare(fs, "Particles");
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "Jets");
      DirectFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      DressedLeptons dressed_leps(fs, bare_leps, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 20*GeV);
      declare(dressed_leps, "Leptons");
      declare(MissingMomentum(fs), "MET");

      // Initialise ROOT file & tree
      _tf = make_unique<TFile>(getOption("ROOTFILE", "Rivet.root").c_str(), "RECREATE");
      _tt = make_unique<TTree>("Rivet", "Rivet physics objects");
      _tt->Branch("npart", &_npart);
      _tt->Branch("nchpart", &_nchpart);
      _tt->Branch("njet", &_njet);
      _tt->Branch("nbjet", &_nbjet);
      _tt->Branch("j1pt", &_j1pt);
      _tt->Branch("j1eta", &_j1eta);
      _tt->Branch("nelec", &_nelec);
      _tt->Branch("e1pt", &_e1pt);
      _tt->Branch("e1eta", &_e1eta);
      _tt->Branch("nmuon", &_nmuon);
      _tt->Branch("m1pt", &_m1pt);
      _tt->Branch("m1eta", &_m1eta);
      _tt->Branch("met", &_met);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve objects
      Particles particles = apply<FinalState>(event, "Particles").particles();
      Particles chparticles = filter_select(particles, isCharged);
      Particles leptons = apply<FinalState>(event, "Leptons").particlesByPt();
      Particles elecs = filter_select(leptons, isElectron);
      Particles muons = filter_select(leptons, isMuon);
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV);
      idiscardIfAnyDeltaRLess(jets, leptons, 0.2);
      Jets bjets = filter_select(jets, hasBTag(Cuts::pT > 5*GeV && Cuts::abseta < 2.5));
      double met = apply<MissingMomentum>(event, "MET").missingPt();

      // Write event row
      _npart = particles.size();
      _nchpart = chparticles.size();
      _njet = jets.size();
      _nbjet = bjets.size();
      _j1pt = (jets.size() > 0) ? jets[0].pT()/GeV : NAN;
      _j1eta = (jets.size() > 0) ? jets[0].eta()/GeV : NAN;
      _nelec = elecs.size();
      _e1pt = (elecs.size() > 0) ? elecs[0].pT()/GeV : NAN;
      _e1eta = (elecs.size() > 0) ? elecs[0].eta()/GeV : NAN;
      _nmuon = muons.size();
      _m1pt = (muons.size() > 0) ? elecs[0].pT()/GeV : NAN;
      _m1eta = (muons.size() > 0) ? elecs[0].eta()/GeV : NAN;
      _met = met;
      _tt->Fill();
    }


    /// Finalize
    void finalize() {
      _tf->Write();
    }

    /// @}


    /// @name Tree variables
    /// @{
    int _npart, _nchpart, _njet, _nbjet, _nelec, _nmuon;
    double _j1pt, _j1eta, _e1pt, _e1eta, _m1pt, _m1eta, _met;
    /// @}


    /// Output file
    unique_ptr<TFile> _tf;
    /// Output tree
    unique_ptr<TTree> _tt;

  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_NTUPLE_ROOT);

}
