// -*- C++ -*-
#include "Rivet/Projections/InvMassFinalState.hh"

namespace Rivet {


  InvMassFinalState::InvMassFinalState(const FinalState& fsp,
                                       const pair<PdgId, PdgId>& idpair, // pair of decay products
                                       double minmass, // min inv mass
                                       double maxmass, // max inv mass
                                       double masstarget)
    : _minmass(minmass), _maxmass(maxmass), _masstarget(masstarget), _useTransverseMass(false)
  {
    setName("InvMassFinalState");
    declare(fsp, "FS");
    _decayids.push_back(idpair);
  }


  InvMassFinalState::InvMassFinalState(const FinalState& fsp,
                                       const vector<pair<PdgId, PdgId> >& idpairs,  // vector of pairs of decay products
                                       double minmass, // min inv mass
                                       double maxmass, // max inv mass
                                       double masstarget)
    : _decayids(idpairs), _minmass(minmass), _maxmass(maxmass), _masstarget(masstarget), _useTransverseMass(false)
  {
    setName("InvMassFinalState");
    declare(fsp, "FS");
  }


  InvMassFinalState::InvMassFinalState(const pair<PdgId, PdgId>& idpair, // pair of decay products
                                       double minmass, // min inv mass
                                       double maxmass, // max inv mass
                                       double masstarget)
    : _minmass(minmass), _maxmass(maxmass), _masstarget(masstarget), _useTransverseMass(false)
  {
    setName("InvMassFinalState");
    _decayids.push_back(idpair);
  }


  InvMassFinalState::InvMassFinalState(const vector<pair<PdgId, PdgId> >& idpairs,  // vector of pairs of decay products
                                       double minmass, // min inv mass
                                       double maxmass, // max inv mass
                                       double masstarget)
    : _decayids(idpairs), _minmass(minmass), _maxmass(maxmass), _masstarget(masstarget), _useTransverseMass(false)
  {
    setName("InvMassFinalState");
  }


  CmpState InvMassFinalState::compare(const Projection& p) const {
    // First compare the final states we are running on
    CmpState fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != CmpState::EQ) return fscmp;

    // Then compare the two as final states
    const InvMassFinalState& other = dynamic_cast<const InvMassFinalState&>(p);
    fscmp = FinalState::compare(other);
    if (fscmp != CmpState::EQ) return fscmp;

    // Compare the mass limits
    CmpState masstypecmp = cmp(_useTransverseMass, other._useTransverseMass);
    if (masstypecmp != CmpState::EQ) return masstypecmp;
    CmpState massllimcmp = cmp(_minmass, other._minmass);
    if (massllimcmp != CmpState::EQ) return massllimcmp;
    CmpState masshlimcmp = cmp(_maxmass, other._maxmass);
    if (masshlimcmp != CmpState::EQ) return masshlimcmp;

    // Compare the decay species
    CmpState decaycmp = cmp(_decayids, other._decayids);
    if (decaycmp != CmpState::EQ) return decaycmp;

    // Finally compare them as final states
    return FinalState::compare(other);
  }


  void InvMassFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    calc(fs.particles());
  }


  void InvMassFinalState::calc(const Particles& inparticles) {
    _theParticles.clear();
    _particlePairs.clear();

    // Containers for the particles of type specified in the pair
    vector<const Particle*> type1, type2;
    // Get all the particles of the type specified in the pair from the particle list
    for (const Particle& ipart : inparticles) {
      // Loop around possible particle pairs
      for (const PdgIdPair& ipair : _decayids) {
        if (ipart.pid() == ipair.first) {
          if (accept(ipart)) type1 += &ipart;
        } else if (ipart.pid() == ipair.second) {
          if (accept(ipart)) type2 += &ipart;
        }
      }
    }
    if (type1.empty() || type2.empty()) return;

    // Temporary container of selected particles iterators
    // Useful to compare iterators and avoid double occurrences of the same
    // particle in case it matches with more than another particle
    vector<const Particle*> tmp;

    // Now calculate the inv mass
    pair<double, pair<Particle, Particle> > closestPair;
    closestPair.first = 1e30;
    for (const Particle* i1 : type1) {
      for (const Particle* i2 : type2) {
        // Check this is actually a pair
        // (if more than one pair in vector particles can be unrelated)
        bool found = false;
        for (const PdgIdPair& ipair : _decayids) {
          if (i1->pid() == ipair.first && i2->pid() == ipair.second) {
            found = true;
            break;
          }
        }
        if (!found) continue;

        FourMomentum v4 = i1->momentum() + i2->momentum();
        if (v4.mass2() < 0) {
          MSG_DEBUG("Constructed negative inv mass2: skipping!");
          continue;
        }
        bool passedMassCut = false;
        if (_useTransverseMass) {
          passedMassCut = inRange(mT(i1->momentum(), i2->momentum()), _minmass, _maxmass);
        } else {
          passedMassCut = inRange(v4.mass(), _minmass, _maxmass);
        }

        if (passedMassCut) {
          MSG_DEBUG("Selecting particles with IDs " << i1->pid() << " & " << i2->pid()
                    << " and mass = " << v4.mass()/GeV << " GeV");
          // Store accepted particles, avoiding duplicates
          if (find(tmp.begin(), tmp.end(), i1) == tmp.end()) {
            tmp.push_back(i1);
            _theParticles += *i1;
          }
          if (find(tmp.begin(), tmp.end(), i2) == tmp.end()) {
            tmp.push_back(i2);
            _theParticles += *i2;
          }
          // Store accepted particle pairs
          _particlePairs += make_pair(*i1, *i2);
          if (_masstarget>0.0) {
            double diff=fabs(v4.mass()-_masstarget);
            if (diff<closestPair.first) {
              closestPair.first=diff;
              closestPair.second=make_pair(*i1, *i2);
            }
          }
        }
      }
    }
    if (_masstarget > 0.0 && closestPair.first < 1e30) {
      _theParticles.clear();
      _particlePairs.clear();
      _theParticles += closestPair.second.first;
      _theParticles += closestPair.second.second;
      _particlePairs += closestPair.second;
    }

    MSG_DEBUG("Selected " << _theParticles.size() << " particles " << "(" << _particlePairs.size() << " pairs)");
    if (getLog().isActive(Log::TRACE)) {
      for (const Particle& p : _theParticles) {
        MSG_TRACE("PDG ID: " << p.pid() << ", ID: " << HepMCUtils::uniqueId(p.genParticle())) ;
      }
    }
  }


  /// Constituent pairs
  const std::vector<std::pair<Particle, Particle> >& InvMassFinalState::particlePairs() const {
    return _particlePairs;
  }


}
