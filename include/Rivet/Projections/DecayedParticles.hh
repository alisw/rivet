// -*- C++ -*-
#ifndef RIVET_DecayedParticles_HH
#define RIVET_DecayedParticles_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Find the decay products of particles in the projection for subsquent analyses
  class DecayedParticles : public Projection {
  public:

    /// @name Standard constructors etc.
    //@{
    /// Constructor.
    DecayedParticles() {}

    DecayedParticles(const ParticleFinder & particles) {
      setName("DecayedParticles");
      declare(particles, "PARTICLES");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(DecayedParticles);
    
    /// Virtual destructor.
    virtual ~DecayedParticles() { }
    //@}

  public :
    
    /// Add a particle to be considered stable when finding the decay products
    DecayedParticles & addStable(PdgId pid) {
      _stable.insert(pid);
      return *this;
    }

    /**
     *  Access to the decaying particles
     */
    const Particles & decaying() const {return _decaying;}

    /**
     *  Access to the number of stable particles
     */
    const vector<unsigned int> & nStable() const {return _nStable;}

    /**
     *  Access to the decay products
     */
    const vector<map<PdgId,Particles> > & decayProducts() const {return _products;}

    /**
     *  Check the particles in the ith mode
     */
    bool modeMatches(size_t imode,unsigned int nstable, map<PdgId,unsigned int> prod) const {
      // same no of stable particles
      if(nstable!=_nStable[imode]) return false;
      for (const auto & kv : prod ) {
	// check if same decay products
	map<PdgId,Particles>::const_iterator iloc = _products[imode].find(kv.first);
	// same type of product
	if (iloc == _products[imode].end()) return false;
	// and same number
	if(iloc->second.size()!=kv.second) return false;
      }
      // pass all the tests
      return true;
    }
    
  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// Compare projections.
    virtual CmpState compare(const Projection& p) const;

  private :

    /**
     * Recursive function to find the decay products
     */ 
    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   map<PdgId,Particles> & products);
    
  private :

    /**
     *  Stable particles
     */
    set<PdgId> _stable;

    /**
     *  The decaying particles
     */
    Particles _decaying;

    /**
     *  The number of stable decay products
     */
    vector<unsigned int> _nStable;

    /**
     *  The decay products
     */
    vector<map<PdgId,Particles> > _products;
  };


}

#endif
