#include "Rivet/Tools/ParticleUtils.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  FirstParticleWith::FirstParticleWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  FirstParticleWithout::FirstParticleWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  LastParticleWith::LastParticleWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  LastParticleWithout::LastParticleWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }


  HasParticleAncestorWith::HasParticleAncestorWith(const Cut& c, bool only_physical)
    : fn([&](const Particle& p){ return c->accept(p); }),
      onlyphysical(only_physical) { }

  HasParticleAncestorWithout::HasParticleAncestorWithout(const Cut& c, bool only_physical)
    : fn([&](const Particle& p){ return c->accept(p); }),
      onlyphysical(only_physical) { }

  HasParticleParentWith::HasParticleParentWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleParentWithout::HasParticleParentWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleChildWith::HasParticleChildWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleChildWithout::HasParticleChildWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleDescendantWith::HasParticleDescendantWith(const Cut& c, bool remove_duplicates)
    : fn([&](const Particle& p){ return c->accept(p); }),
      rmduplicates(remove_duplicates) { }

  HasParticleDescendantWithout::HasParticleDescendantWithout(const Cut& c, bool remove_duplicates)
    : fn([&](const Particle& p){ return c->accept(p); }),
      rmduplicates(remove_duplicates) { }


  Particles& ifilter_select(Particles& particles, const Cut& c) {
    if (c == Cuts::OPEN) return particles;
    // return ifilter_select(particles, *c);
    return ifilter_select(particles, [&](const Particle& p){return c->accept(p);});
  }

  Particles& ifilter_discard(Particles& particles, const Cut& c) {
    if (c == Cuts::OPEN) { particles.clear(); return particles; }
    // return ifilter_discard(particles, *c);
    return ifilter_discard(particles, [&](const Particle& p){return c->accept(p);});
  }


  /// @brief Check if the list of particles given is compatible with the requested list of pids
  /// @note if abspid is true, then only absolute values of the Particles' pids are compared.
  bool partsAre(const Particles& originalparts, const vector<int>& originalpids, bool abspid) {

    // obviously if the inputs have different lengths, they aren't compatible.
    if (originalparts.size() != originalpids.size())
      return false;

    vector<int> pids = originalpids;
    sort(pids.begin(), pids.end());

    Particles parts = originalparts;

    // sort the Particles by abspid or pid, as desired.
    if (abspid) {
      sort
      ( parts.begin()
      , parts.end()
      , [](const Particle& p1, const Particle& p2) { return p1.abspid() < p2.abspid(); }
      );
    }
    else {
      sort
      ( parts.begin()
      , parts.end()
      , [](const Particle& p1, const Particle& p2) { return p1.pid() < p2.pid(); }
      );
    }

    // once sorted, every pid should align exactly;
    // otherwise these aren't the parts we're looking for!
    for (size_t i = 0; i < pids.size(); i++) {
      if (abspid) {
        if (parts[i].abspid() != pids[i]) return false;
      } else {
        if (parts[i].pid() != pids[i]) return false;
      }
    }

    return true;
  }

  // checks if the decay chains of the input Particles are compatible with the
  // desired pids.
  // this takes a list of particles rather than a single particle so that it can
  // be called recursively.
  // this algorithm probably scales pretty badly with the parts and pids
  // lengths.
  /// @note initial photons are _not_ ignored!
  /// @todo this is really some kind of tree structure, but in this implementation there are duplicate "branches" that just cause inefficiency... but it shouldn't be _wrong_.
  
  bool cascadeContains(const Particles& parts, const vector<int>& pids, bool absolute, bool ignorephoton) {

    if (parts.size() > pids.size()) {
      // base case 1: can't reduce the number of particles to the target size.
      return false;

    } else if (parts.size() == pids.size()) {
        // base case 2: this is the right number of particles. Are they the ones
        // we're looking for?
        return partsAre(parts, pids, absolute);

    } else {
      // recursive case: parts.size() < pids.size(), which means there's still a
      // chance that they match further along in the decay chain. 
      // for each particle left in our list, we expand it to its decay products
      // and see if it matches what we're looking for (possibly with subsequent expansions).
      for (size_t i = 0; i < parts.size(); i++) {
        Particles originalchildren = parts[i].children();

        // empty children -> start over
        if (!originalchildren.size())
          continue;

        // build the list of viable children of particle i, ignoring photons if
        // requested.
        Particles children;
        if (ignorephoton) {
          for (size_t j = 0; j < originalchildren.size(); j++) {
            const Particle& part = originalchildren[j];

            if (part.pid() != PID::PHOTON) {
              children.push_back(part);
            }

          }
        } else
          children = originalchildren;

        // build the new "current" set of particles, with the decay of particle
        // i expanded to its children.
        const Particles newparts =
          slice(parts, 0, i) + children + slice(parts, i+1, parts.size());

        // recursive call: check if the "current" set of particles fits our
        // needs (possibly after further decays).
        if (cascadeContains(newparts, pids, absolute, ignorephoton))
          return true;
      }

    }

    return false;
  }

}
