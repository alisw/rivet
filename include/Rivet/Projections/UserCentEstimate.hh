// -*- C++ -*-
#ifndef RIVET_USERCENTESTIMATE_HH
#define RIVET_USERCENTESTIMATE_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"


namespace Rivet {

class UserCentEstimate: public SingleValueProjection {
public:

  UserCentEstimate() {
    setName("UserCentEstimate");
    declare(HepMCHeavyIon(), "HepMC");
  }

  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(UserCentEstimate);

protected:

  void project(const Event& e) {
    clear();
    set(apply<HepMCHeavyIon>(e, "HepMC").user_cent_estimate());
   }

  CmpState compare(const Projection& p) const {
    return CmpState::EQ;
  }

};

}

#endif
