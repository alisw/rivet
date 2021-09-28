// -*- C++ -*-
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {


/// Generic analysis looking at various distributions of final state particles
class ALICE_2015_PPBCentrality : public Analysis {

public:

  DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2015_PPBCentrality);

  /// Book histograms and initialise projections before the run
  void init() {

    // One projection for the actual observable, and one for the
    // generated impact parameter.
    declare(ALICE::V0AMultiplicity(), "V0A");
    declare(ImpactParameterProjection(), "IMP");
    declare(ALICE::V0AndTrigger(), "Trigger");

    // The calibration histogram:
    book(_calib, "V0A", 100, 0.0, 500.0);


    // The alternative histogram based on impact parameter. Note that
    // it MUST be named the same as the histogram for the experimental
    // observable with an added _IMP suffix.
    book(_impcalib, "V0A_IMP", 400, 0.0, 20.0);

  }

  /// Perform the per-event analysis
  void analyze(const Event& event) {

    // The alternative centrality based on generated impact
    // parameter, assumes that the generator does not describe the
    // full final state, and should therefore be filled even if the
    // event is not triggered.
    _impcalib->fill(apply<SingleValueProjection>(event, "IMP")());

    if ( !apply<ALICE::V0AndTrigger>(event, "Trigger")() ) vetoEvent;

    _calib->fill(apply<ALICE::V0AMultiplicity>(event, "V0A")());

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
DECLARE_RIVET_PLUGIN(ALICE_2015_PPBCentrality);

}
