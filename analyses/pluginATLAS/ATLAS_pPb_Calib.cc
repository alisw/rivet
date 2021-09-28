// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/AtlasCommon.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"

namespace Rivet {


/// Generic analysis looking at various distributions of final state particles
class ATLAS_pPb_Calib : public Analysis {

public:

  DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_pPb_Calib);

  /// Book histograms and initialise projections before the run
  void init() {

    // One projection for the actual observable, and one for the
    // generated impact parameter.
    declare(ATLAS::SumET_PB_Centrality(), "Centrality");
    declare(ImpactParameterProjection(), "IMP");
    declare(ATLAS::MinBiasTrigger(), "Trigger");

    // The calibrationhistogram:
    book(_calib, "SumETPb", 100, 0.0, 200.0);

    // The alternative histogram based on impact parameter. Note that
    // it MUST be named the same as the histogram for the experimental
    // observable with an added _IMP suffix for the Pecentile<>
    // binning to work properly.
    book(_impcalib, "SumETPb_IMP", 400, 0.0, 20.0);


  }
  
  /// Perform the per-event analysis
  void analyze(const Event& event) {
    
    // The alternative centrality based on generated impact
    // parameter, assumes that the generator does not describe the
    // full final state, and should therefore be filled even if the
    // event is not triggered.
    _impcalib->fill(apply<SingleValueProjection>(event, "IMP")());

    if ( !apply<ATLAS::MinBiasTrigger>(event, "Trigger")() ) vetoEvent;

    _calib->fill(apply<ATLAS::SumET_PB_Centrality>(event, "Centrality")());

  }
  
  /// Finalize
  void finalize() {

    _calib->normalize();
    _impcalib->normalize();

  }

private:

  /// The calibration histograms.
  Histo1DPtr _calib;
  Histo1DPtr _impcalib;

};


// The hook for the plugin system
DECLARE_RIVET_PLUGIN(ATLAS_pPb_Calib);

}
