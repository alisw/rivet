// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class L3_1990_I298078 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_1990_I298078);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections to use
      const FinalState FS;
      declare(FS, "FS");
      FastJets jadeJets = FastJets(FS, FastJets::JADE, 0.7, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jadeJets, "JadeJets");
      // book histos
      book(_h_y_2_JADE, 1,1,1);
      book(_h_y_3_JADE, 1,1,2);
      book(_h_y_4_JADE, 1,1,3);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FastJets& jadejet = apply<FastJets>(event, "JadeJets");
      if (jadejet.clusterSeq()) {
	const double y_23 = jadejet.clusterSeq()->exclusive_ymerge_max(2);
	const double y_34 = jadejet.clusterSeq()->exclusive_ymerge_max(3);
	const double y_45 = jadejet.clusterSeq()->exclusive_ymerge_max(4);
	for (size_t i = 0; i < _h_y_2_JADE->numBins(); ++i) {
	  double ycut = _h_y_2_JADE->bin(i).xMid();
	  double width = _h_y_2_JADE->bin(i).width();
	  if (y_23 < ycut) _h_y_2_JADE->fillBin(i,width);
	}
	for (size_t i = 0; i < _h_y_3_JADE->numBins(); ++i) {
	  double ycut = _h_y_3_JADE->bin(i).xMid();
	  double width = _h_y_3_JADE->bin(i).width();
	  if (y_34 < ycut && y_23 > ycut) {
	    _h_y_3_JADE->fillBin(i,width);
	  }
	}
	for (size_t i = 0; i < _h_y_4_JADE->numBins(); ++i) {
	  double ycut = _h_y_4_JADE->bin(i).xMid();
	  double width = _h_y_4_JADE->bin(i).width();
	  if (y_45 < ycut && y_34 > ycut) {
	    _h_y_4_JADE->fillBin(i,width);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_y_2_JADE, 100./ sumOfWeights());
      scale(_h_y_3_JADE, 100./ sumOfWeights());
      scale(_h_y_4_JADE, 100./ sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_y_2_JADE,_h_y_3_JADE,_h_y_4_JADE;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(L3_1990_I298078);

}
