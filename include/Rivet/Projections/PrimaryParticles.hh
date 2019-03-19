// -*- C++ -*-
/**
 * @file   PrimaryParticles.hh
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Mon Aug 27 08:46:34 2018
 * @brief  Primary particle definition based on PDG code. 
 */
#ifndef RIVET_PrimaryParticles_HH
#define RIVET_PrimaryParticles_HH

#include <Rivet/Projections/ParticleFinder.hh>
#include <Rivet/Tools/Cuts.hh>

namespace Rivet
{
  /// @brief Project out primary particles according to definition.
  /// A Rivet projection that mimics an experimental primary partcile
  /// definition by projecting out according to particle id.
  /// The projection can be further specialized to accomodate 
  /// specific experimental definitions.
  //
  class PrimaryParticles : public ParticleFinder
  {
  public:
    /** 
     * Constructor 
     * 
     * @param cuts   Normal particle cuts 
     * @param pdgIds List of PDG IDs which are considered primary. 
     *
     * @todo Instead of using a real vector use an initializer list -
     * more flexible.  Also, do not provide a default for the PDG IDs
     * (even if empty) as this will force the user to actively select
     * which codes are primary.
     */
    PrimaryParticles(std::initializer_list<int> pdgIds,
		     const Cut& c=Cuts::open()) :
      ParticleFinder(c), _pdgIds(pdgIds) {
      setName("PrimaryParticles");
    }
    // Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PrimaryParticles);
      
    /**
     * Copy constructor    
     *
     * @param other Other object to copy from 
     */
    PrimaryParticles(const PrimaryParticles& other) : 
      ParticleFinder(other), _pdgIds(other._pdgIds) {
    }
    /** 
     * Compare to projections.  
     * 
     * @param p Projection to compare to. 
     * 
     * @return true (EQUIVALENT) if the projection @a p is of the same
     * type as this, if the cuts are equal, and that the list of PDG
     * IDs are the same.
     */
    virtual int compare(const Projection& p) const
    {
      const PrimaryParticles* other = dynamic_cast<const PrimaryParticles*>(&p);
      if (!other) return UNDEFINED;
      if (_cuts != other->_cuts || _pdgIds != other->_pdgIds) return UNDEFINED;
      return EQUIVALENT;
      
    }
  protected:
    /**
     * Do the projection.
     *
     * @param e Event to project from 
     */
    virtual void project(const Event& e);

    /** 
     * @{
     * @name Internally used member functions 
     */
    /** 
     * Check if the particle is a priamry. 
     *
     * @param p Pointer to a HepMC particle 
     *
     * @return true if the particle @a p is considered primary 
     */
    virtual bool isPrimary(const HepMC::GenParticle* p) const;
    /**
     * Check if the particle should be ignored by the status code of
     * the particle.
     */
    virtual bool isIgnored(const HepMC::GenParticle* p) const;
    /** 
     * Check PDG ID of particle @a p is in the list of accepted
     * primaries.
     *
     * @param p Particle to investigate.
     *
     * @return true if the particle PDG ID is in the list of known
     * primary PDG IDs.
     */
    virtual bool isPrimaryPID(const HepMC::GenParticle* p) const;
    /*
     * Check if a particle @a p has decayed.
     *
     * @param p Pointer to HepMC particle 
     *
     * @return true if the particle has decayed according to the
     * status flag of the particle @a p
     */
    virtual bool hasDecayed(const HepMC::GenParticle* p) const;
    /**
     * Check if a particle is a beam (remnant) particle.
     *
     * @param p Particle to check 
     *
     * @return true if the particle @a p is a (remnant) beam particle 
     */
    virtual bool isBeam(const HepMC::GenParticle* p) const;
    /*
     * Get the immediate ancestor of a particle.
     * 
     * @param p Particle for which to get the immediate ancestor 
     *
     * @return Pointer to immediate ancestor or null if there's no ancestor. 
     */
    const HepMC::GenParticle* ancestor(const HepMC::GenParticle* p) const;
    /*
     * Get the immediate ancestor of a particle, which is @e not an
     * ignored particle.
     * 
     * @param p Particle for which to get the immediate ancestor 
     *
     * @return Pointer to immediate ancestor or null if there's no ancestor. 
     */
    const HepMC::GenParticle* ancestor(const HepMC::GenParticle* p, bool) const;
    /** Particle types to test for */
    std::vector<int> _pdgIds;
  };
}
#endif
