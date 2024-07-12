// -*- C++ -*-
#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  CmpState InitialQuarks::compare(const Projection& p) const {
    return CmpState::EQ;
  }


  void InitialQuarks::project(const Event& e) {
    _theParticles.clear();

    for (ConstGenParticlePtr p : HepMCUtils::particles(e.genEvent())) {
      ConstGenVertexPtr pv = p->production_vertex();
      ConstGenVertexPtr dv = p->end_vertex();
      const PdgId pid = abs(p->pdg_id());
      bool passed = inRange((long)pid, 1, 6);
      if (passed) {
        if (pv != 0) {
          for (ConstGenParticlePtr pp : HepMCUtils::particles(pv, Relatives::PARENTS)){
            // Only accept if parent is electron or Z0
            const PdgId pid = abs(pp->pdg_id());
            passed = (pid == PID::ELECTRON || abs(pp->pdg_id()) == PID::ZBOSON || abs(pp->pdg_id()) == PID::GAMMA);
          }
        } else {
          passed = false;
        }
      }

      if (getLog().isActive(Log::TRACE)) {
        const int st = p->status();
        const double pT = p->momentum().perp();
        const double eta = p->momentum().eta();
        MSG_TRACE(std::boolalpha
                  << "ID = " << p->pdg_id() << ", status = " << st << ", pT = " << pT
                  << ", eta = " << eta << ": result = " << passed);
        if (pv != 0) {
          for (ConstGenParticlePtr pp : HepMCUtils::particles(pv, Relatives::PARENTS)) {
            MSG_TRACE(std::boolalpha << " parent ID = " << pp->pdg_id());
          }
        }
        if (dv != 0) {
          for (ConstGenParticlePtr pp : HepMCUtils::particles(dv, Relatives::CHILDREN)) {
            MSG_TRACE(std::boolalpha << " child ID  = " << pp->pdg_id());
          }
        }
      }
      if (passed) _theParticles.push_back(Particle(p));
    }
    MSG_DEBUG("Number of initial quarks = " << _theParticles.size());
    if (!_theParticles.empty()) {
      for (size_t i = 0; i < _theParticles.size(); i++) {
        MSG_DEBUG("Initial quark[" << i << "] = " << _theParticles[i].pid());
      }
    }
  }


}
