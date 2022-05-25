// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief jet rates at 91 GeV
  class DELPHI_1990_I297698 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1990_I297698);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(FastJets(cfs, FastJets::JADE, 0.7), "JadeJets");
      // histos
      book(_h_2,1,1,1);
      book(_h_3,1,1,2);
      book(_h_4,1,1,3);
      book(_h_5,1,1,4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      const FastJets& jets = apply<FastJets>(event, "JadeJets");
      if (jets.clusterSeq()) {
        const double y_23 = jets.clusterSeq()->exclusive_ymerge_max(2);
        const double y_34 = jets.clusterSeq()->exclusive_ymerge_max(3);
        const double y_45 = jets.clusterSeq()->exclusive_ymerge_max(4);
        const double y_56 = jets.clusterSeq()->exclusive_ymerge_max(5);
        for (size_t i = 0; i < _h_2->numBins(); ++i) {
          double ycut = _h_2->bin(i).xMid();
          double width = _h_2->bin(i).xWidth();
          if (y_23 < ycut) {
            _h_2->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_3->numBins(); ++i) {
          double ycut = _h_3->bin(i).xMid();
          double width = _h_3->bin(i).xWidth();
          if (y_34 < ycut && y_23 > ycut) {
            _h_3->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_4->numBins(); ++i) {
          double ycut = _h_4->bin(i).xMid();
          double width = _h_4->bin(i).xWidth();
          if (y_45 < ycut && y_34 > ycut) {
            _h_4->fillBin(i, width);
          }
        }
        for (size_t i = 0; i < _h_5->numBins(); ++i) {
          double ycut = _h_5->bin(i).xMid();
          double width = _h_5->bin(i).xWidth();
          if (y_56 < ycut && y_45 > ycut) {
            _h_5->fillBin(i, width);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_2, 100/sumOfWeights());
      scale(_h_3, 100/sumOfWeights());
      scale(_h_4, 100/sumOfWeights());
      scale(_h_5, 100/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_2,_h_3,_h_4,_h_5;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(DELPHI_1990_I297698);

}
