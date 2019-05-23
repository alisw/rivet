// -*- C++ -*-
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  void UnstableParticles::project(const Event& e) {
    _theParticles.clear();

    /// @todo Replace PID veto list with PID:: functions?
    vector<PdgId> vetoIds;
    vetoIds += 22; // status 2 photons don't count!
    vetoIds += 110; vetoIds += 990; vetoIds += 9990; // Reggeons
    for (const GenParticle* p : Rivet::particles(e.genEvent())) {
      const Particle rp(p);
      const int st = p->status();
      bool passed =
        (st == 1 || (st == 2 && !contains(vetoIds, abs(p->pdg_id())))) &&
        !PID::isParton(p->pdg_id()) && ///< Always veto partons
        !p->is_beam() && // Filter beam particles
        _cuts->accept(rp);

      // Avoid double counting by re-marking as unpassed if ID == (any) parent ID
      const GenVertex* pv = p->production_vertex();
      // if (passed && pv) {
      //   for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin(); pp != pv->particles_in_const_end(); ++pp) {
      //     if (p->pdg_id() == (*pp)->pdg_id() && (*pp)->status() == 2) {
      //       passed = false;
      //       break;
      //     }
      //   }
      // }
      //

      // Avoid double counting by re-marking as unpassed if ID == any child ID
      const GenVertex* dv = p->end_vertex();
      if (passed && dv) {
        for (GenParticle* pp : particles_out(dv)) {
          if (p->pdg_id() == pp->pdg_id() && pp->status() == 2) {
            passed = false;
            break;
          }
        }
      }

      // Add to output particles collection
      if (passed) _theParticles.push_back(rp);

      // Log parents and children
      if (getLog().isActive(Log::TRACE)) {
        MSG_TRACE("ID = " << p->pdg_id()
                  << ", status = " << st
                  << ", pT = " << p->momentum().perp()
                  << ", eta = " << p->momentum().eta()
                  << ": result = " << std::boolalpha << passed);
        if (pv) {
          for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin(); pp != pv->particles_in_const_end(); ++pp) {
            MSG_TRACE("  parent ID = " << (*pp)->pdg_id());
          }
        }
        if (dv) {
          for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
            MSG_TRACE("  child ID  = " << (*pp)->pdg_id());
          }
        }
      }
    }
    MSG_DEBUG("Number of unstable final-state particles = " << _theParticles.size());
  }


}
