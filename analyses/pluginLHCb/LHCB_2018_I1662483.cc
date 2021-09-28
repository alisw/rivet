// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Measurement of forward top pair production in the mu e channel at LHCb
  ///
  /// Both leptons with pt > 20 GeV and in pseudorapidity range 2.0 - 4.5.
  /// Leptons are separated from each other by a dR of 0.1.
  /// Leading jet is required to be a b-jet with pt > 20 GeV and pseudorapidity in range 2.2 - 4.2.
  /// The jet is separated from both leptons by dR of 0.5.
  ///
  class LHCB_2018_I1662483 : public Analysis {
  public:

    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LHCB_2018_I1662483);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection:
      // all final-state particles
      const FinalState fs;

      //get the final state b quarks for b-tagging
      FinalPartons bquarks(Cuts::abspid == 5 && Cuts::pT > 5*GeV);
      declare(bquarks, "bquarks");

      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.5
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.5, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // FinalState of bare muons and electrons in the event
      Cut lepton_cuts = Cuts::pT>=20*GeV && Cuts::etaIn(2.0, 4.5);
      FinalState fsmuons(Cuts::abspid == PID::MUON && lepton_cuts );
      FinalState fselectrons(Cuts::abspid == PID::ELECTRON && lepton_cuts);

      declare(fsmuons, "muons");
      declare(fselectrons, "electrons");

      // Inclusive cross-section measurement
      book(_h_fiducial_xsect ,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // retrieve leptons, sorted by pT
      const FinalState& muons = apply<FinalState>(event, "muons");
      const FinalState& electrons = apply<FinalState>(event, "electrons");
      const FinalPartons& bquarks = apply<FinalPartons>(event, "bquarks");

      // retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::etaIn(2.2,4.2));

      if (jets.empty()) vetoEvent;
      bool pass = false;
      for (const auto& m : muons.particles()) {
        for (const auto& e : electrons.particles()) {
          if (deltaR(m.momentum(), e.momentum()) < 0.1 ) continue;
          vector<Jet> lepton_jets;
          for (const auto& j : jets ){
            if (deltaR(j.momentum(), m.momentum()) > 0.5 &&
                deltaR(j.momentum(), e.momentum()) > 0.5)
              lepton_jets.push_back(j);
            if (lepton_jets.size() > 0){
              for (const auto& b : bquarks.particles()) {
                if ( deltaR(b.momentum(), lepton_jets.at(0).momentum()) < 0.5 )  {
                  pass = true;
                }
              }
            }
          }
        }
      }

      // veto event if it doesn't pass our selection
      if (!pass) vetoEvent;
      // append to cross-section
      _h_fiducial_xsect->fill(sqrtS()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_fiducial_xsect, crossSection()/femtobarn/sumOfWeights()); // norm to cross section
    }


    /// Histogram
    Histo1DPtr _h_fiducial_xsect;

  };


  DECLARE_RIVET_PLUGIN(LHCB_2018_I1662483);

}
