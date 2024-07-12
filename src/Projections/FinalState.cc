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
      doubles offShellM2s, vtxDisps;
      for (ConstGenParticlePtr p : HepMCUtils::particles(e.genEvent())) {
        if (p->status() == 1) {
          MSG_TRACE("FS GV = " << p->production_vertex()->position());
          Particle rp(p);
          /// @todo Complete off-shell testing with comparison to a dict of pole masses
          const double m2 = rp.mass2();
          if (m2 < -1*GeV2) offShellM2s += m2; //< only note significantly negative m2s for now
          /// @todo Reinstate vertex displacement warnings with a more refined primary-particles
          /// definition: this screams about ctau0 > 10mm SM particles, which is not helpful
          // const double disp = rp.origin().polarRadius();
          // if (disp > 10*cm) vtxDisps += disp;
          // if (disp > 100*cm) MSG_WARNING("Particle origin transverse-displaced by more than 1m: " << rp << " at " << rp.origin());
          _theParticles.push_back(rp);
        }
      }
      MSG_TRACE("Number of open-FS selected particles = " << _theParticles.size());
      if (!offShellM2s.empty()) MSG_WARNING(offShellM2s.size() << " final-state particles found with negative mass^2:" << offShellM2s);
      if (!vtxDisps.empty()) MSG_WARNING(vtxDisps.size() << " final-state particles found with significantly transverse-displaced origin vertices:" << vtxDisps);
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
