// -*- C++ -*-
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  int MissingMomentum::compare(const Projection& p) const {
    return mkNamedPCmp(p, "VisibleFS");
  }


  void MissingMomentum::clear() {
    _momentum = FourMomentum();
    _set = 0.0;
    _vet = Vector3();
  }


  void MissingMomentum::project(const Event& e) {
    clear();

    // Project into final state
    const FinalState& vfs = applyProjection<FinalState>(e, "VisibleFS");
    foreach (const Particle& p, vfs.particles()) {
      const FourMomentum& mom = p.momentum();
      _momentum += mom;
      _set += mom.Et();
      _vet += mom.Et() * mom.vector3().setZ(0.0).unit();
    }
  }


  const FourMomentum MissingMomentum::visibleMomentum(double mass) const {
    /// @todo Couldn't we just reset the internal _momentum's mass and return by value? Would require mutable, though
    FourMomentum p4 = _momentum;
    const double pmod2 = p4.p3().mod2();
    const double new_energy = sqrt(pmod2 + sqr(mass));
    p4.setE(new_energy);
    return p4;
  }


}
