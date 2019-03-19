#include "Rivet/Particle.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"
#include "Rivet/Tools/ParticleUtils.hh"

namespace Rivet {


  void Particle::setConstituents(const Particles& cs, bool setmom) {
    _constituents = cs;
    if (setmom) _momentum = sum(cs, p4, FourMomentum());
  }


  void Particle::addConstituent(const Particle& c, bool addmom) {
    _constituents += c;
    if (addmom) _momentum += c;
  }


  void Particle::addConstituents(const Particles& cs, bool addmom) {
    _constituents += cs;
    if (addmom)
      for (const Particle& c : cs)
        _momentum += c;
  }


  Particles Particle::rawConstituents() const {
    if (!isComposite()) return Particles{*this};
    Particles rtn;
    for (const Particle& p : constituents()) rtn += p.rawConstituents();
    return rtn;
  }


  Particle& Particle::transformBy(const LorentzTransform& lt) {
    _momentum = lt.transform(_momentum);
    return *this;
  }


  bool Particle::isVisible() const {
    // Charged particles are visible
    if ( PID::threeCharge(pid()) != 0 ) return true;
    // Neutral hadrons are visible
    if ( PID::isHadron(pid()) ) return true;
    // Photons are visible
    if ( pid() == PID::PHOTON ) return true;
    // Gluons are visible (for parton level analyses)
    if ( pid() == PID::GLUON ) return true;
    // Everything else is invisible
    return false;
  }


  bool Particle::isStable() const {
    return genParticle() != NULL &&
      genParticle()->status() == 1 &&
      genParticle()->end_vertex() == NULL;
  }


  Particles Particle::ancestors(const Cut& c, bool physical_only) const {
    Particles rtn;

    // this case needed protecting against (at least for the latest Herwig... not sure why
    // it didn't show up earlier
    if (genParticle() == NULL) return rtn;

    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->production_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::children))
    //   rtn += Particle(gp);
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::ancestors); it != gv->particles_end(HepMC::ancestors); ++it) {
      if (physical_only && (*it)->status() != 1 && (*it)->status() != 2) continue;
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  Particles Particle::parents(const Cut& c) const {
    Particles rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->production_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::children))
    //   rtn += Particle(gp);
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::parents); it != gv->particles_end(HepMC::parents); ++it) {
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  Particles Particle::children(const Cut& c) const {
    Particles rtn;
    /// @todo Something going wrong with taus -> GenParticle nullptr?
    if (genParticle() == nullptr) return rtn;
    if (isStable()) return rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    //cout << genParticle()->end_vertex() << endl;
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->end_vertex() );
    if (gv == nullptr) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::children))
    //   rtn += Particle(gp);
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::children); it != gv->particles_end(HepMC::children); ++it) {
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
  /// @todo Use recursion through replica-avoiding functions to avoid bookkeeping duplicates
  Particles Particle::allDescendants(const Cut& c, bool remove_duplicates) const {
    Particles rtn;
    if (isStable()) return rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->end_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::descendants))
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it != gv->particles_end(HepMC::descendants); ++it) {
      const Particle p(*it);
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      if (remove_duplicates && (*it)->end_vertex() != NULL) {
        // size_t n = 0; ///< @todo Only remove 1-to-1 duplicates?
        bool dup = false;
        /// @todo Yuck, HepMC
        for (GenVertex::particle_iterator it2 = (*it)->end_vertex()->particles_begin(HepMC::children); it2 != (*it)->end_vertex()->particles_end(HepMC::children); ++it2) {
          // n += 1; if (n > 1) break;
          if ((*it)->pdg_id() == (*it2)->pdg_id()) { dup = true; break; }
        }
        if (dup) continue;
      }
      rtn += p;
    }
    return rtn;
  }


  /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
  Particles Particle::stableDescendants(const Cut& c) const {
    Particles rtn;
    if (isStable()) return rtn;
    /// @todo Remove this const mess crap when HepMC doesn't suck
    GenVertexPtr gv = const_cast<GenVertexPtr>( genParticle()->end_vertex() );
    if (gv == NULL) return rtn;
    /// @todo Would like to do this, but the range objects are broken
    // foreach (const GenParticlePtr gp, gv->particles(HepMC::descendants))
    for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it != gv->particles_end(HepMC::descendants); ++it) {
      // if ((*it)->status() != 1 || (*it)->end_vertex() != NULL) continue;
      const Particle p(*it);
      if (!p.isStable()) continue;
      if (c != Cuts::OPEN && !c->accept(p)) continue;
      rtn += p;
    }
    return rtn;
  }


  double Particle::flightLength() const {
    if (isStable()) return -1;
    if (genParticle() == NULL) return 0;
    if (genParticle()->production_vertex() == NULL) return 0;
    const HepMC::FourVector v1 = genParticle()->production_vertex()->position();
    const HepMC::FourVector v2 = genParticle()->end_vertex()->position();
    return sqrt(sqr(v2.x()-v1.x()) + sqr(v2.y()-v1.y()) + sqr(v2.z()-v1.z()));
  }


  bool Particle::hasParent(PdgId pid) const {
    return hasParentWith(hasPID(pid));
  }

  bool Particle::hasParentWith(const Cut& c) const {
    return hasParentWith([&](const Particle& p){return c->accept(p);});
  }


  bool Particle::hasAncestor(PdgId pid, bool only_physical) const {
    return hasAncestorWith(hasPID(pid), only_physical);
  }

  bool Particle::hasAncestorWith(const Cut& c, bool only_physical) const {
    return hasAncestorWith([&](const Particle& p){return c->accept(p);}, only_physical);
  }


  bool Particle::hasChildWith(const Cut& c) const {
    return hasChildWith([&](const Particle& p){return c->accept(p);});
  }


  bool Particle::hasDescendantWith(const Cut& c, bool remove_duplicates) const {
    return hasDescendantWith([&](const Particle& p){return c->accept(p);}, remove_duplicates);
  }

  bool Particle::hasStableDescendantWith(const Cut& c) const {
    return hasStableDescendantWith([&](const Particle& p){return c->accept(p);});
  }



  bool Particle::fromBottom() const {
    return hasAncestorWith([](const Particle& p){
        return p.genParticle()->status() == 2 && p.isHadron() && p.hasBottom();
      });
  }

  bool Particle::fromCharm() const {
    return hasAncestorWith([](const Particle& p){
        return p.genParticle()->status() == 2 && p.isHadron() && p.hasCharm();
      });
  }

  bool Particle::fromHadron() const {
    return hasAncestorWith([](const Particle& p){
        return p.genParticle()->status() == 2 && p.isHadron();
      });
  }

  bool Particle::fromTau(bool prompt_taus_only) const {
    if (prompt_taus_only && fromHadron()) return false;
    return hasAncestorWith([](const Particle& p){
        return p.genParticle()->status() == 2 && isTau(p);
      });
  }

  bool Particle::fromHadronicTau(bool prompt_taus_only) const {
    return hasAncestorWith([&](const Particle& p){
        return p.genParticle()->status() == 2 && isTau(p) && (!prompt_taus_only || p.isPrompt()) && hasHadronicDecay(p);
      });
  }


  bool Particle::isDirect(bool allow_from_direct_tau, bool allow_from_direct_mu) const {
    while (!_isDirect.second) { ///< @todo Replace awkward caching with C++17 std::optional
      // Immediate short-circuit: hadrons can't be direct, and for partons we can't tell
      if (isHadron() || isParton()) {
        _isDirect = std::make_pair(false, true); break;
      }

      // Obtain links to parentage
      if (genParticle() == nullptr) { _isDirect = std::make_pair(false, true); break; } // no HepMC connection, give up! Throw UserError exception?
      const GenVertexPtr prodVtx = genParticle()->production_vertex();
      if (prodVtx == nullptr) { _isDirect = std::make_pair(false, true); break; } // orphaned particle, has to be assume false
      const pair<GenParticlePtr, GenParticlePtr> beams = prodVtx->parent_event()->beam_particles();

      /// @todo Would be nicer to be able to write this recursively up the chain, exiting as soon as a parton or string/cluster is seen
      for (const GenParticlePtr ancestor : Rivet::particles(prodVtx, HepMC::ancestors)) {
        const PdgId pid = ancestor->pdg_id();
        if (ancestor->status() != 2) continue; // no non-standard statuses or beams to be used in decision making
        if (ancestor == beams.first || ancestor == beams.second) continue; // PYTHIA6 uses status 2 for beams, I think... (sigh)
        if (PID::isParton(pid)) continue; // PYTHIA6 also uses status 2 for some partons, I think... (sigh)
        if (PID::isHadron(pid)) { _isDirect = std::make_pair(false, true); break; } // direct particles can't be from hadron decays
        if (abs(pid) == PID::TAU && abspid() != PID::TAU && !allow_from_direct_tau) { _isDirect = std::make_pair(false, true); break; } // allow or ban particles from tau decays (permitting tau copies)
        if (abs(pid) == PID::MUON && abspid() != PID::MUON && !allow_from_direct_mu) { _isDirect = std::make_pair(false, true); break; } // allow or ban particles from muon decays (permitting muon copies)
      }
      if (!_isDirect.second) _isDirect = std::make_pair(true, true); //< guarantee loop exit
    }
    return _isDirect.first;
  }



  ///////////////////////


  // DISABLED UNTIL VANILLA CC7 COMPATIBILITY NOT NEEDED

  // /// Particles copy constructor from vector<Particle>
  // Particles::Particles(const std::vector<Particle>& vps) : base(vps) {}

  // /// Particles -> FourMomenta cast/conversion operator
  // Particles::operator FourMomenta () const {
  //   // FourMomenta rtn(this->begin(), this->end());
  //   FourMomenta rtn; rtn.reserve(this->size());
  //   for (size_t i = 0; i < this->size(); ++i) rtn.push_back((*this)[i]);
  //   return rtn;
  // }

  // /// Particles concatenation operator
  // Particles operator + (const Particles& a, const Particles& b) {
  //   Particles rtn(a);
  //   rtn += b;
  //   return rtn;
  // }


  //////////////////////////////////


  /// Allow a Particle to be passed to an ostream.
  std::ostream& operator << (std::ostream& os, const Particle& p) {
    string pname;
    try {
      pname = PID::toParticleName(p.pid());
    } catch (...) {
      pname = "PID=" + to_str(p.pid());
    }
    os << "Particle<" << pname << " @ " << p.momentum()/GeV << " GeV>";
    return os;
  }


  /// Allow ParticlePair to be passed to an ostream.
  std::ostream& operator << (std::ostream& os, const ParticlePair& pp) {
    os << "[" << pp.first << ", " << pp.second << "]";
    return os;
  }



}
