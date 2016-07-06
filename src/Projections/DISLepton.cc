// -*- C++ -*-
#include "Rivet/Projections/DISLepton.hh"

namespace Rivet {


  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = pcast<DISLepton>(p);
    return mkNamedPCmp(other, "Beam") || mkNamedPCmp(other, "FS");
  }


  void DISLepton::project(const Event& e) {
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();

    Particle inLepton;

    bool firstIsLepton = PID::isLepton(inc.first.pid());
    bool secondIsLepton = PID::isLepton(inc.second.pid());

    if (firstIsLepton && !secondIsLepton) {
      _incoming = inc.first;
    } else if (!firstIsLepton && secondIsLepton) {
      _incoming = inc.second;
    } else {
      //eek!
      throw	Error("DISLepton projector could not find the correct beam.");
    }

    _sign = (_incoming.momentum().pz() > 0.0)? 1.0: -1.0;
    long id = _incoming.pid();

    double pzMax = -1000000000.0;

    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    foreach (const Particle& p, fs.particles()) {
      double pz = _sign * p.momentum().pz();
      if (p.pid() == id && pz > pzMax) {
        _outgoing = p;
        pzMax = pz;
      }
    }

    if (_outgoing.genParticle() == NULL) {
      throw Error("DISLepton projector could not find the scattered lepton.");
    }
  }


}
