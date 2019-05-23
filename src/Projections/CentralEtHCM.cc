// -*- C++ -*-
#include "Rivet/Projections/CentralEtHCM.hh"

namespace Rivet {


  void CentralEtHCM::project(const Event& e) {
    const DISFinalState& fs = applyProjection<DISFinalState>(e, "FS");
    if ( fs.failed() ) {
      fail();
      return;
    }
    _sumet = 0.0;
    for (const Particle& p : fs.particles()) {
      /// @todo Generalise rapidity cut value
      if (fabs(p.rapidity()) < 0.5) _sumet += p.Et();
    }
  }


}
