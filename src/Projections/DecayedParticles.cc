// -*- C++ -*-
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {
  
  CmpState DecayedParticles::compare(const Projection& p) const {
    const DecayedParticles& other = dynamic_cast<const DecayedParticles&>(p);
    // Compare particles definitions
    const CmpState teq = mkPCmp(other, "PARTICLES");
    if (teq != CmpState::EQ) return teq;
    // Compare set of stable particles 
    const CmpState nfeq = cmp(_stable.size(), other._stable.size());
    if (nfeq != CmpState::EQ) return nfeq;
    for(const PdgId & pid : _stable) {
      if (other._stable.find(pid)==other._stable.end()) return CmpState::NEQ;
    }
    // If we got this far, we're equal
    return CmpState::EQ;
  }
  
  void DecayedParticles::findDecayProducts(const Particle & mother, unsigned int & nstable,
					   map<PdgId,Particles> & products) {
    for(const Particle & p : mother.children()) {
      PdgId pid = p.pid();
      // no decay products or in list of stable particles
      if ( p.children().empty() || _stable.find(pid)!=_stable.end()) {
	++nstable;
	map<PdgId,Particles>::iterator iloc = products.find(pid);
	if(iloc!=products.end()) iloc->second.push_back(p);
	else products[pid] = Particles({p});
      }	
      else
	findDecayProducts(p,nstable,products);
    }
  }


  /// Perform the particle and decay product finding
  void DecayedParticles::project(const Event& e) {
    Particles part = apply<ParticleFinder>(e, "PARTICLES").particles();
    // clear the storage and get the new particles
    _decaying.clear(); _decaying.reserve(part.size());
    _nStable.clear();  _nStable.reserve(part.size());
    _products.clear(); _products.reserve(part.size());
    for (const Particle& p : part) {
      if(p.children().size()<=1) continue;
      _decaying.push_back(p);
      unsigned int nstable(0);
      map<PdgId,Particles> products;
      findDecayProducts(p,nstable,products);
      _nStable.push_back(nstable);
      _products.push_back(products);
    }
  }
  
}
