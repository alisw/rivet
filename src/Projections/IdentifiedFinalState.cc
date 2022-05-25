// -*- C++ -*-
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  IdentifiedFinalState::IdentifiedFinalState(const FinalState& fsp, const vector<PdgId>& pids) {
    setName("IdentifiedFinalState");
    declare(fsp, "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const FinalState& fsp, PdgId pid) {
    setName("IdentifiedFinalState");
    declare(fsp, "FS");
    acceptId(pid);
  }

  IdentifiedFinalState::IdentifiedFinalState(const Cut& c, const vector<PdgId>& pids) {
    setName("IdentifiedFinalState");
    declare(FinalState(c), "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const vector<PdgId>& pids, const Cut& c) {
    setName("IdentifiedFinalState");
    declare(FinalState(c), "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const Cut& c, PdgId pid) {
    setName("IdentifiedFinalState");
    declare(FinalState(c), "FS");
    acceptId(pid);
  }

  IdentifiedFinalState::IdentifiedFinalState(PdgId pid, const Cut& c) {
    setName("IdentifiedFinalState");
    declare(FinalState(c), "FS");
    acceptId(pid);
  }



  CmpState IdentifiedFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != CmpState::EQ) return fscmp;

    const IdentifiedFinalState& other = dynamic_cast<const IdentifiedFinalState&>(p);
    CmpState pidssize = cmp(_pids.size(), other._pids.size());
    if (pidssize != CmpState::EQ) return pidssize;
    return cmp(_pids, other._pids);
  }


  void IdentifiedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    _remainingParticles.clear();
    _remainingParticles.reserve(fs.particles().size());
    for (const Particle& p : fs.particles()) {
      if (acceptedIds().find(p.pid()) != acceptedIds().end()) {
        _theParticles.push_back(p);       // Identified
      } else {
        _remainingParticles.push_back(p); // Remaining
      }
    }
  }


}
