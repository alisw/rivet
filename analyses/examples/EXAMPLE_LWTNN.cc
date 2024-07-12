// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetLWTNN.hh"

namespace Rivet {


  /// @brief Example analysis demonstrating use of LWTNN neural networks
  ///
  /// An example analysis to demonstrate using LWTNN to preserve a neural net
  /// inside Rivet. The network used here is a "fake" model that takes in
  /// info about the 3 leading-pT jets and returns an arbitrary score that
  /// is plotted.
  class EXAMPLE_LWTNN : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_LWTNN);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetfs, "jets");

      // Initialise the LWTNN graph.
      // Doing this in init saves loading it separately for each event.
      _lwg = mkLWTNN(analysisDataPath("json"));

      // Book a histograms with custom binning
      book(_h["DNN_output"], "DNN_output", 20, 0.0, 100.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);

      // Require at least three jets:
      if (jets.size() < 3){
        vetoEvent;
      }

      // Fill your map of variables for input to the lwtnn DNN
      map<string, double> inputs = {{"jet1_pt", jets[0].pT()},
                           {"jet2_pt", jets[1].pT()},
                           {"jet3_pt", jets[2].pT()},
                           {"jet1_eta", jets[0].eta()},
                           {"jet2_eta", jets[1].eta()},
                           {"jet3_eta", jets[2].eta()}};
      map<string, double> output = _lwg->compute(inputs);

      // Fill histogram with leading b-jet pT
      _h["DNN_output"]->fill(output["output"]);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Nothing to do here -- unique_ptr does our lwtnn clean up for us.
    }

    /// @}


    /// Histograms
    map<string, Histo1DPtr> _h;

    /// LWTNN neural-network pointer
    unique_ptr<lwt::LightweightNeuralNetwork> _lwg;

  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_LWTNN);

}
