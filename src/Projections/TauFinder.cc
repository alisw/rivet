// -*- C++ -*-
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  void TauFinder::project(const Event& e) {
    _theParticles.clear();
    const auto& ufs = applyProjection<UnstableParticles>(e, "UFS");
    for (const Particle& p : ufs.particles()) {
      if (p.abspid() != PID::TAU) continue;
      if (_dectype == ANY || (_dectype == LEPTONIC && isLeptonic(p)) || (_dectype == HADRONIC && isHadronic(p)) )
        _theParticles.push_back(p);
    }
  }


  int TauFinder::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "UFS");
    if (fscmp != EQUIVALENT) return fscmp;

    const auto& other = dynamic_cast<const TauFinder&>(p);
    return cmp(_dectype, other._dectype);
  }


}
