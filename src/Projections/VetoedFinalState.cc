// -*- C++ -*-
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  int VetoedFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;
    /// @todo We can do better than this...
    if (_vetofsnames.size() != 0) return UNDEFINED;
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    return \
      cmp(_vetoCuts, other._vetoCuts) ||
      cmp(_compositeVetoes, other._compositeVetoes) ||
      cmp(_parentVetoes, other._parentVetoes);
  }


  void VetoedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());

    // Veto by cut
    if (_vetoCuts.empty()) {
      _theParticles = fs.particles();
    } else {
      // Test every particle against the codes
      for (const Particle& p : fs.particles()) {
        bool match = false;
        for (const Cut& c : _vetoCuts) {
          if (c->accept(p)) {
            match = true;
            break;
          }
        }
        if (!match) {
          MSG_TRACE("Storing particle " << p);
          _theParticles.push_back(p);
        } else {
          MSG_TRACE("Vetoing particle " << p);
        }
      }
    }

    // Composite particle mass vetoing
    /// @todo YUCK! Clean up...
    set<Particles::iterator> toErase;
    for (auto nIt = _nCompositeDecays.begin(); nIt != _nCompositeDecays.end() && !_theParticles.empty(); ++nIt) {
      map<set<Particles::iterator>, FourMomentum> oldMasses;
      map<set<Particles::iterator>, FourMomentum> newMasses;
      set<Particles::iterator> start;
      start.insert(_theParticles.begin());
      oldMasses.insert(pair<set<Particles::iterator>, FourMomentum>(start, _theParticles.begin()->momentum()));
      for (int nParts = 1; nParts != *nIt; ++nParts) {
        for (auto mIt = oldMasses.begin(); mIt != oldMasses.end(); ++mIt) {
          Particles::iterator pStart = *(mIt->first.rbegin());
          for (auto pIt = pStart + 1; pIt != _theParticles.end(); ++pIt) {
            FourMomentum cMom = mIt->second + pIt->momentum();
            set<Particles::iterator> pList(mIt->first);
            pList.insert(pIt);
            newMasses[pList] = cMom;
          }
        }
        oldMasses = newMasses;
        newMasses.clear();
      }
      for (auto mIt = oldMasses.begin(); mIt != oldMasses.end(); ++mIt) {
        double mass2 = mIt->second.mass2();
        if (mass2 >= 0.0) {
          double mass = sqrt(mass2);
          for (auto cIt = _compositeVetoes.lower_bound(*nIt); cIt != _compositeVetoes.upper_bound(*nIt); ++cIt) {
            auto massRange = cIt->second;
            if (mass < massRange.second && mass > massRange.first) {
              for (auto lIt = mIt->first.begin(); lIt != mIt->first.end(); ++lIt) {
                toErase.insert(*lIt);
              }
            }
          }
        }
      }
    }
    for (auto p = toErase.rbegin(); p != toErase.rend(); ++p) {
      _theParticles.erase(*p);
    }

    // Remove particles whose parents match entries in the parent veto PDG ID codes list
    /// @todo There must be a nice way to do this -- an STL algorithm (or we provide a nicer wrapper)
    for (PdgId vetoid : _parentVetoes) {
      for (Particles::iterator ip = _theParticles.begin(); ip != _theParticles.end(); ++ip) {
        const GenVertex* startVtx = ip->genParticle()->production_vertex();
        if (startVtx == NULL) continue;
        // Loop over parents and test their IDs
        /// @todo Could use any() here?
        for (const GenParticle* parent : Rivet::particles(startVtx, HepMC::ancestors)) {
          if (vetoid == parent->pdg_id()) {
            ip = _theParticles.erase(ip); --ip; //< Erase this _theParticles entry
            break;
          }
        }
      }
    }

    // Finally veto on the registered FSes
    for (const string& ifs : _vetofsnames) {
      const ParticleFinder& vfs = applyProjection<ParticleFinder>(e, ifs);
      const Particles& pvetos = vfs.rawParticles();
      ifilter_discard(_theParticles, [&](const Particle& pcheck) {
          if (pcheck.genParticle() == nullptr) return false;
          for (const Particle& pveto : pvetos) {
            if (pveto.genParticle() == nullptr) continue;
            if (pveto.genParticle() == pcheck.genParticle()) { MSG_TRACE("Vetoing: " << pcheck); return true; }
          }
          return false;
        });
    }

    MSG_DEBUG("FS vetoing from #particles = " << fs.size() << " -> " << _theParticles.size());
  }


}
