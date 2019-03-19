// -*- C++ -*-
#ifndef RIVET_USERCENTESTIMATE_HH
#define RIVET_USERCENTESTIMATE_HH

#include "Rivet/Projections/SingleValueProjection.hh"


namespace Rivet {

class UserCentEstimate: public SingleValueProjection {

public:
  
  UserCentEstimate() {
    setName("UserCentEstimate");
  }

  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(UserCentEstimate);

protected:

  void project(const Event& e) {
    clear();
#if HEPMC_VERSION_CODE >= 3000000
    const HepMC::HeavyIon * hi = e.genEvent()->heavy_ion();
    if ( hi && hi->user_cent_estimate >= 0.0 )
      set(hi->centrality*100.0);
#endif
   }
  
  int compare(const Projection& p) const {
    return 0;
  }
  
};

}

#endif

