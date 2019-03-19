// -*- C++ -*-
#ifndef RIVET_IMPACTPARAMETERPROJECTION_HH
#define RIVET_IMPACTPARAMETERPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"


namespace Rivet {

class ImpactParameterProjection: public SingleValueProjection {

public:
  
  ImpactParameterProjection() {
    setName("ImpactParameterProjection");
  }

  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(ImpactParameterProjection);

protected:

  void project(const Event& e) {
    clear();
    const HepMC::HeavyIon * hi = e.genEvent()->heavy_ion();
    if ( hi && hi->impact_parameter() >= 0.0 )
      set(hi->impact_parameter());
  }
  
  int compare(const Projection& p) const {
    return 0;
  }
  
};

}

#endif

