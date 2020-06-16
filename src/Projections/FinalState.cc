// -*- C++ -*-
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  FinalState::FinalState(const Cut& c)
    : ParticleFinder(c)
  {
    setName("FinalState");
    const bool isopen = (c == Cuts::open());
    MSG_TRACE("Check for open FS conditions: " << std::boolalpha << isopen);
    if (!isopen) declare(FinalState(), "OpenFS");
  }


  FinalState::FinalState(const FinalState& fsp, const Cut& c)
    : ParticleFinder(c)
  {
    setName("FinalState");
    MSG_TRACE("Registering base FSP as 'PrevFS'");
    declare(fsp, "PrevFS");
  }


  CmpState FinalState::compare(const Projection& p) const {
    const FinalState& other = dynamic_cast<const FinalState&>(p);
    // First check if there is a PrevFS and it it matches
    if (hasProjection("PrevFS") != other.hasProjection("PrevFS")) return CmpState::NEQ;
    if (hasProjection("PrevFS")) {
      const PCmp prevcmp = mkPCmp(other, "PrevFS");
      if (prevcmp != CmpState::EQ) return  CmpState::NEQ;
    }
    // Then check the extra cuts
    const bool cutcmp = _cuts == other._cuts;
    MSG_TRACE(_cuts << " VS " << other._cuts << " -> EQ == " << std::boolalpha << cutcmp);
    if (!cutcmp) return CmpState::NEQ;
    // Checks all passed: these FSes are equivalent
    return CmpState::EQ;
  }


  void FinalState::project(const Event& e) {
    _theParticles.clear();

    // Handle "open FS" special case, which should not/cannot recurse
    if (_cuts == Cuts::OPEN) {
      MSG_TRACE("Open FS processing: should only see this once per event (" << e.genEvent()->event_number() << ")");
      for (ConstGenParticlePtr p : HepMCUtils::particles(e.genEvent())) {
        if (p->status() == 1) {
          MSG_TRACE("FS GV = " << p->production_vertex()->position());
          _theParticles.push_back(Particle(p));
        }
      }
      MSG_TRACE("Number of open-FS selected particles = " << _theParticles.size());
      return;
    }

    // Base the calculation on PrevFS if available, otherwise OpenFS
    /// @todo In general, we'd like to calculate a restrictive FS based on the most restricted superset FS.
    const Particles& allstable = applyProjection<FinalState>(e, (hasProjection("PrevFS") ? "PrevFS" : "OpenFS")).particles();
    MSG_TRACE("Beginning Cuts selection");
    for (const Particle& p : allstable) {
      const bool passed = accept(p);
      MSG_TRACE("Choosing: ID = " << p.pid()
                << ", pT = " << p.pT()/GeV << " GeV"
                << ", eta = " << p.eta()
                << ": result = " << std::boolalpha << passed);
      if (passed) _theParticles.push_back(p);
    }
    MSG_TRACE("Number of final-state particles = " << _theParticles.size());
  }


  /// Decide if a particle is to be accepted or not.
  bool FinalState::accept(const Particle& p) const {
    // Not having status == 1 should never happen!
    assert(p.genParticle() == NULL || p.genParticle()->status() == 1);
    return _cuts->accept(p);
  }


}
