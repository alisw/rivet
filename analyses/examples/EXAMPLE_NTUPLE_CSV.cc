// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include <fstream>

namespace Rivet {


  /// @brief An example of writing a CSV text file with columnar data per event
  class EXAMPLE_NTUPLE_CSV : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_NTUPLE_CSV);


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

      // Initialise CSV file
      _fout = std::ofstream(getOption("CSVFILE", "Rivet.csv"));
      writefield("npart", 5);
      writefield("nchpart", 7);
      writefield("njet", 5);
      writefield("nbjet", 5);
      writefield("j1pt", 5);
      writefield("j1eta", 5);
      writefield("nelec", 5);
      writefield("e1pt", 5);
      writefield("e1eta", 5);
      writefield("nmuon", 5);
      writefield("m1pt", 5);
      writefield("m1eta", 5);
      writefield("met", 5, "");
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

      // Write CSV row
      writefield(particles.size(), 5);
      writefield(chparticles.size(), 7);
      writefield(jets.size(), 5);
      writefield(bjets.size(), 5);
      writefield(jets.size() > 0 ? jets[0].pT()/GeV : NAN, 5);
      writefield(jets.size() > 0 ? jets[0].eta()/GeV : NAN, 5);
      writefield(elecs.size(), 5);
      writefield(elecs.size() > 0 ? elecs[0].pT()/GeV : NAN, 5);
      writefield(elecs.size() > 0 ? elecs[0].eta()/GeV : NAN, 5);
      writefield(muons.size(), 5);
      writefield(muons.size() > 0 ? elecs[0].pT()/GeV : NAN, 5);
      writefield(muons.size() > 0 ? elecs[0].eta()/GeV : NAN, 5);
      writefield(met, 5, "");
    }


    /// Finalize (nothing to do since ofstream closes automatically)
    // void finalize() {}

    /// @}


    /// Helper method for writing columns in the CSV
    template <typename T>
    void writefield(const T& x, size_t width, const string& post=", ") {
      _fout.precision(width);
      _fout << std::setw(width) << x << post;
      if (post.empty()) _fout << endl;
    }


    /// CSV output file
    std::ofstream _fout;

  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_NTUPLE_CSV);

}
