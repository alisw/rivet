// -*- C++ -*-
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  FinalState::FinalState(const Cut& c)
    : ParticleFinder(c)
  {
    setName("FinalState");
    const bool isopen = (c == Cuts::open());
    MSG_TRACE("Check for open FS conditions: " << std::boolalpha << isopen);
    if (!isopen) addProjection(FinalState(), "OpenFS");
  }


  FinalState::FinalState(const FinalState& fsp, const Cut& c)
    : ParticleFinder(c)
  {
    setName("FinalState");
    MSG_TRACE("Registering base FSP as 'PrevFS'");
    addProjection(fsp, "PrevFS");
  }


  /// @deprecated, keep for backwards compatibility for now.
  FinalState::FinalState(double mineta, double maxeta, double minpt) {
    setName("FinalState");
    const bool openpt = isZero(minpt);
    const bool openeta = (mineta <= -MAXDOUBLE && maxeta >= MAXDOUBLE);
    MSG_TRACE("Check for open FS conditions:" << std::boolalpha << " eta=" << openeta << ", pt=" << openpt);
    if (openpt && openeta) {
      _cuts = Cuts::open();
    } else {
      addProjection(FinalState(), "OpenFS");
      if (openeta)
        _cuts = (Cuts::pT >= minpt);
      else if ( openpt )
        _cuts = Cuts::etaIn(mineta, maxeta);
      else
        _cuts = (Cuts::etaIn(mineta, maxeta) && Cuts::pT >= minpt);
    }
  }


  int FinalState::compare(const Projection& p) const {
    const FinalState& other = dynamic_cast<const FinalState&>(p);
    // First check if there is a PrevFS and it it matches
    if (hasProjection("PrevFS") != other.hasProjection("PrevFS")) return UNDEFINED;
    if (hasProjection("PrevFS")) {
      const PCmp prevcmp = mkPCmp(other, "PrevFS");
      if (prevcmp != EQUIVALENT) return prevcmp;
    }
    // Then check the extra cuts
    const bool cutcmp = _cuts == other._cuts;
    MSG_TRACE(_cuts << " VS " << other._cuts << " -> EQ == " << std::boolalpha << cutcmp);
    if (!cutcmp) return UNDEFINED;
    // Checks all passed: these FSes are equivalent
    return EQUIVALENT;
  }


  void FinalState::project(const Event& e) {
    _theParticles.clear();

    // Handle "open FS" special case, which should not/cannot recurse
    if (_cuts == Cuts::OPEN) {
      MSG_TRACE("Open FS processing: should only see this once per event (" << e.genEvent()->event_number() << ")");
      for (const GenParticle* p : Rivet::particles(e.genEvent())) {
        if (p->status() == 1) {
          MSG_TRACE("FS GV = " << p->production_vertex());
          _theParticles.push_back(Particle(*p));
        }
      }
      return;
    }

    // Base the calculation on PrevFS if available, otherwise OpenFS
    /// @todo In general, we'd like to calculate a restrictive FS based on the most restricted superset FS.
    const Particles& allstable = applyProjection<FinalState>(e, (hasProjection("PrevFS") ? "PrevFS" : "OpenFS")).particles();
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
