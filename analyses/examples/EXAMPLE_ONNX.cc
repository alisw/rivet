// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetONNXrt.hh"

namespace Rivet {


  /// @brief Example analysis demonstrating use of ONNX neural networks
  ///
  /// An example analysis to demonstrate using ONNX to preserve a neural net
  /// inside Rivet. The network used here is a "fake" model that takes in
  /// info about the 3 leading-pT jets and returns an arbitrary score that
  /// is plotted.
  /// This is very similar to EXAMPLE_LWTNN but the network is not the same.
  /// (though once again, the output of the network is completely meaningless)
  class EXAMPLE_ONNX : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_ONNX);

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

      // Initialise our neural net.
      // Much quicker to load in init than once per event!
      _nn = make_unique<RivetONNXrt>(RivetONNXrt(analysisDataPath("onnx")));
      // ^This finds a file with the same name as the analysis but with the
      // 'onnx' suffix (i.e. EXAMPLE_ONNX.onnx).

      // Book a histogram with custom binning
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

      // Create vector of inputs.
      vector<float> nn_input = {
        (float)jets[0].pt(),
        (float)jets[1].pt(),
        (float)jets[2].pt(),
        (float)jets[0].eta(),
        (float)jets[1].eta(),
        (float)jets[2].eta()
      };
      // Vector to fill with output
      vector<float> nn_output;

      // Compute
      _nn->compute(nn_input, nn_output);

      // Fill histogram
      _h["DNN_output"]->fill(abs(nn_output[0]));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Nothing to do here.
    }

    /// @}

    
    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    /// @}

    private:
    /// @name Member variables
    /// @{
    /// The neural network
    unique_ptr<RivetONNXrt> _nn;
    /// @}
  };


  RIVET_DECLARE_PLUGIN(EXAMPLE_ONNX);
}
