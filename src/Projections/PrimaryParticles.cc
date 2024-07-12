#include <Rivet/Projections/PrimaryParticles.hh>
#include <Rivet/Particle.hh>
#include <Rivet/Event.hh>
#include <Rivet/Tools/ParticleIdUtils.hh>
#include <Rivet/Tools/ParticleName.hh>

namespace Rivet {


  void PrimaryParticles::project(const Event& e)
  {
    _theParticles.clear(); // Clear cache
    bool open = _cuts == Cuts::open(); 
    for (auto p : HepMCUtils::particles(e.genEvent())) {
      if (isPrimary(p) && (open || _cuts->accept(Particle(p)))) 
	_theParticles.push_back(Particle(*p));
    }
  }

  bool PrimaryParticles::isPrimary(ConstGenParticlePtr p) const
  {
    if(isIgnored(p))  return false;
    if(!isPrimaryPID(p)) return false;

    // Loop back over ancestors that are _not_ ignored 
    ConstGenParticlePtr m = p;
    while ((m = ancestor(m,true))) {
      if (isBeam(m)) return true;
      if (isPrimaryPID(m)) return false;
      if (!hasDecayed(m)) return false;
    }
    return true;
  }
  bool PrimaryParticles::isIgnored(ConstGenParticlePtr p) const
  {
    return (p->status()==0 || (p->status()>=11 && p->status()<=200));
  }
    
  bool PrimaryParticles::isPrimaryPID(ConstGenParticlePtr p) const
  {
    int thisPID = abs(p->pdg_id());
    for(const auto pid : _pdgIds)
      if (thisPID == pid) return true;
    return false;
  }

  bool PrimaryParticles::isBeam(ConstGenParticlePtr p) const
  {
    // Pythia6 uses 3 for initial state
    return p && (p->status() == 3 || p->status() == 4); 
  }

  bool PrimaryParticles::hasDecayed(ConstGenParticlePtr p) const
  {
    return p && p->status() == 2;
  }

  ConstGenParticlePtr
  PrimaryParticles::ancestor(ConstGenParticlePtr p) const
  {
    ConstGenVertexPtr vtx = p->production_vertex();
    if (!vtx) return 0;
    
    auto parents = HepMCUtils::particles(vtx, Relatives::PARENTS);
    if ( parents.empty() ) return nullptr;
    return parents[0];
  }
  
  ConstGenParticlePtr
  PrimaryParticles::ancestor(ConstGenParticlePtr p, bool) const
  {
    ConstGenParticlePtr m = p;
    do {
      m = ancestor(m);
    } while (m && isIgnored(m));
    return m;
  }
}
