// -*- C++ -*-
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  ChargedFinalState::ChargedFinalState(const FinalState& fsp) {
    setName("ChargedFinalState");
    declare(fsp, "FS");
  }

  ChargedFinalState::ChargedFinalState(const Cut& c) {
    setName("ChargedFinalState");
    declare(FinalState(c), "FS");
  }

  CmpState ChargedFinalState::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }
}

/*namespace {
  inline bool chargedParticleFilter(const Rivet::Particle& p) {
    return Rivet::PID::charge3(p.pid()) == 0;
  }
}*/

namespace Rivet {
  void ChargedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    // _theParticles.clear();
    // std::remove_copy_if(fs.particles().begin(), fs.particles().end(),
    //                     std::back_inserter(_theParticles),
    //                     [](const Rivet::Particle& p) { return p.charge3() == 0; });
    _theParticles = filter_select(fs.particles(), isCharged);
    MSG_DEBUG("Number of charged final-state particles = " << _theParticles.size());
    if (getLog().isActive(Log::TRACE)) {
      for (const Particle&  p : _theParticles) {
        MSG_TRACE("Selected: " << p.pid() << ", charge = " << p.charge());
      }
    }
  }


}
