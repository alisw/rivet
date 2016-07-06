#include "Rivet/Particle.hh"
#include "Rivet/Tools/RivetBoost.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {


  /// @todo Neaten this up with C++11, via one walker function and several uses with lamba tests


  bool Particle::hasAncestor(PdgId pdg_id) const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      if (ancestor->pdg_id() == pdg_id) return true;
    }
    return false;
  }


  bool Particle::fromBottom() const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() == 2 && (PID::isHadron(pid) && PID::hasBottom(pid))) return true;
    }
    return false;
  }


  bool Particle::fromCharm() const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() == 2 && (PID::isHadron(pid) && PID::hasCharm(pid) && !PID::hasBottom(pid))) return true;
    }
    return false;
  }



  bool Particle::fromHadron() const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() == 2 && PID::isHadron(pid)) return true;
    }
    return false;
  }


  bool Particle::fromTau(bool prompt_taus_only) const {
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    if (prompt_taus_only && fromHadron()) return false;
    GenVertex* prodVtx = genParticle()->production_vertex();
    if (prodVtx == NULL) return false;
    foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() == 2 && abs(pid) == PID::TAU) return true;
    }
    return false;
  }


  // bool Particle::fromDecay() const {
  //   /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
  //   GenVertex* prodVtx = genParticle()->production_vertex();
  //   if (prodVtx == NULL) return false;
  //   foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
  //     const PdgId pid = ancestor->pdg_id();
  //     if (ancestor->status() == 2 && (PID::isHadron(pid) || abs(pid) == PID::TAU)) return true;
  //   }
  //   return false;
  // }


}
