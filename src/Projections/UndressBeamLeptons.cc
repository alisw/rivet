// -*- C++ -*-
#include "Rivet/Projections/UndressBeamLeptons.hh"

namespace Rivet {


  CmpState UndressBeamLeptons::compare(const Projection & p) const {
    const UndressBeamLeptons & o =
      dynamic_cast<const UndressBeamLeptons &>(p);
    return cmp(_thetamax, o._thetamax) || mkNamedPCmp(o, "FS");
  }


  void UndressBeamLeptons::project(const Event& e) {
    Beam::project(e);
    if ( _thetamax <= 0.0 ) return;
    bool l1 = _theBeams.first.isChargedLepton();
    bool l2 = _theBeams.second.isChargedLepton();
    if ( !l1 && !l2 ) return;
    FourMomentum b1 = _theBeams.first.momentum();
    FourMomentum b2 = _theBeams.second.momentum();
    Vector3 b10 = b1.vector3();
    Vector3 b20 = b2.vector3();
    for ( auto p : apply<FinalState>(e, "FS").particles() ) {
      if ( p.pid() != PID::PHOTON ) continue;
      if ( p.momentum().angle(b10) < _thetamax )
        b1 -= p.momentum();
      if ( p.momentum().angle(b20) < _thetamax )
        b2 -= p.momentum();
    }
    _theBeams.first.setMomentum(b1);
    _theBeams.second.setMomentum(b2);
  }


}
