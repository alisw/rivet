/**
 * @file   PrimaryParticles.cc
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Mon Aug 27 09:06:03 2018
 * 
 * @brief  Primary particle definition based on PDG IDs. 
 */

#include <Rivet/Projections/PrimaryParticles.hh>
#include <Rivet/Particle.hh>
#include <Rivet/Event.hh>
#include <Rivet/Tools/ParticleIdUtils.hh>
#include <Rivet/Tools/ParticleName.hh>
#include <HepMC/GenParticle.h>
#include <HepMC/GenVertex.h>

namespace Rivet {

  void PrimaryParticles::project(const Event& e)
  {
    _theParticles.clear(); // Clear cache
    bool open = _cuts == Cuts::open(); 
    for (auto p : Rivet::particles(e.genEvent())) {
      if (isPrimary(p) && (open || _cuts->accept(Particle(p)))) 
	_theParticles.push_back(Particle(*p));
    }
  }

  bool PrimaryParticles::isPrimary(const HepMC::GenParticle* p) const
  {
    if(isIgnored(p))  return false;
    if(!isPrimaryPID(p)) return false;

    // Loop back over ancestors that are _not_ ignored 
    const HepMC::GenParticle* m = p;
    while ((m = ancestor(m,true))) {
      if (isBeam(m)) return true;
      if (isPrimaryPID(m)) return false;
      if (!hasDecayed(m)) return false;
    }
    return true;
  }
  bool PrimaryParticles::isIgnored(const HepMC::GenParticle* p) const
  {
    return (p->status()==0 || (p->status()>=11 && p->status()<=200));
  }
    
  bool PrimaryParticles::isPrimaryPID(const HepMC::GenParticle* p) const
  {
    int thisPID = PID::abspid(p->pdg_id());
    for(const auto pid : _pdgIds)
      if (thisPID == pid) return true;
    return false;
  }

  bool PrimaryParticles::isBeam(const HepMC::GenParticle* p) const
  {
    // Pythia6 uses 3 for initial state
    return p && (p->status() == 3 || p->status() == 4); 
  }

  bool PrimaryParticles::hasDecayed(const HepMC::GenParticle* p) const
  {
    return p && p->status() == 2;
  }

  const HepMC::GenParticle*
  PrimaryParticles::ancestor(const HepMC::GenParticle* p) const
  {
    const HepMC::GenVertex* vtx = p->production_vertex();
    if (!vtx) return 0;
    
    HepMC::GenVertex::particles_in_const_iterator i =
      vtx->particles_in_const_begin();
    if (i == vtx->particles_in_const_end()) return 0;
    return *i;
  }
  
  const HepMC::GenParticle*
  PrimaryParticles::ancestor(const HepMC::GenParticle* p, bool) const
  {
    const HepMC::GenParticle* m = p;
    do {
      m = ancestor(m);
    } while (m && isIgnored(m));
    return m;
  }
}
//
// EOF
//
