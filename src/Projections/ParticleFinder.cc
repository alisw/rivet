// -*- C++ -*-
#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {

  /// @todo HOW DO WE COMPARE CUTS OBJECTS?
  CmpState ParticleFinder::compare(const Projection& p) const {
    const ParticleFinder& other = dynamic_cast<const ParticleFinder&>(p);
    return _cuts == other._cuts ? CmpState::EQ : CmpState::NEQ;
  }

}
