// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"

namespace Rivet {


  class ALICE_2015_PPCentrality : public Analysis {
  public:
     DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2015_PPCentrality);

     // Book histograms, initialize projections.
    void init() {
      declare(ALICE::V0AndTrigger(),"V0-AND");
      declare(ALICE::V0MMultiplicity(),"V0M");
      declare(ImpactParameterProjection(), "IMP");
      book(_v0m, "V0M",100,0,200);
      book(_imp, "V0M_IMP",100,0,20);
    }


    // Per-event analysis
    void analyze(const Event& event) {
      _imp->fill(apply<SingleValueProjection>(event,"IMP")());

      // Check if we have any hit in either V0-A or -C.  If not, the
      // event is not selected and we get out.
      if (!apply<ALICE::V0AndTrigger>(event,"V0-AND")()) return;

      // Fill in the V0 multiplicity for this event
      _v0m->fill(apply<ALICE::V0MMultiplicity>(event,"V0M")());
    }


    void finalize() {
      _v0m->normalize();
      _imp->normalize();
    }


    Histo1DPtr _v0m;
    Histo1DPtr _imp;
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2015_PPCentrality);

}
