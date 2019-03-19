// -*- C++ -*-
#ifndef RIVET_GENERATEDPERCENTILEPROJECTION_HH
#define RIVET_GENERATEDPERCENTILEPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"


namespace Rivet {

class GeneratedPercentileProjection: public SingleValueProjection {

public:
  
  GeneratedPercentileProjection() {
    setName("GeneratedPercentileProjection");
  }

  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(GeneratedPercentileProjection);

protected:

  void project(const Event& e) {
    clear();
#if HEPMC_VERSION_CODE >= 3000000
    const HepMC::HeavyIon * hi = e.genEvent()->heavy_ion();
    if ( hi && hi->centrality >= 0.0 )
      set(hi->centrality*100.0);
#endif
   }
  
  int compare(const Projection& p) const {
    return 0;
  }
  
};

}

#endif

