// -*- C++ -*-
#ifndef RIVET_PrimaryParticles_HH
#define RIVET_PrimaryParticles_HH

#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  /// @brief Project out primary particles according to definition
  ///
  /// A Rivet projection that mimics an experimental primary partcile
  /// definition by projecting out according to particle ID.
  /// The projection can be further specialized to accomodate
  /// specific experimental definitions.
  ///
  /// @author Christian Holm Christensen <cholm@nbi.dk>
  class PrimaryParticles : public ParticleFinder {
  public:

    /// Constructor
	///
	/// @param pids  List of PDG IDs which are considered primary
	/// @param c     Normal particle cuts
    PrimaryParticles(std::initializer_list<int> pids,
		     const Cut& c=Cuts::open()) :
      ParticleFinder(c), _pdgIds(pids) {
      setName("PrimaryParticles");
    }

    // Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PrimaryParticles);

    /// Copy constructor
    PrimaryParticles(const PrimaryParticles& other) :
      ParticleFinder(other), _pdgIds(other._pdgIds) {
    }

    /// Compare to another projection
	///
	/// @param p Projection to compare to.
	///
	/// @return Equivalent if the projection @a p is of the same type as this,
	/// the cuts are equal, and the list of PDG IDs is the same.
    virtual CmpState compare(const Projection& p) const
    {
      const PrimaryParticles* other = dynamic_cast<const PrimaryParticles*>(&p);
      if (!other) return CmpState::NEQ;
      if (_cuts != other->_cuts || _pdgIds != other->_pdgIds) return CmpState::NEQ;
      return CmpState::EQ;
    }


  protected:

    /// Do the projection.
	///
	/// @param e Event to project from
    virtual void project(const Event& e);

	/// @name Internally used member functions
    /// @{
    ///
    /// Check if the particle is a primary.
	///
	/// @param p Pointer to a HepMC particle
	///
	/// @return true if the particle @a p is considered primary
    virtual bool isPrimary(ConstGenParticlePtr p) const;

    /// Check if the particle should be ignored, via its status code
    virtual bool isIgnored(ConstGenParticlePtr p) const;

    /// Check PDG ID of particle @a p is in the list of accepted primaries
	///
	/// @param p Particle to investigate.
	///
	/// @return true if the particle PDG ID is in the list of known primary PDG IDs.
    virtual bool isPrimaryPID(ConstGenParticlePtr p) const;

    /// Check if a particle @a p has decayed.
	///
	/// @param p Pointer to HepMC particle
	///
	/// @return true if the particle has decayed according to the
	/// status flag of the particle @a p
    virtual bool hasDecayed(ConstGenParticlePtr p) const;

    /// Check if a particle is a beam (remnant) particle.
	///
	/// @param p Particle to check
	///
	/// @return true if the particle @a p is a (remnant) beam particle
    virtual bool isBeam(ConstGenParticlePtr p) const;

	/// Get the immediate ancestor of a particle.
	///
	/// @param p Particle for which to get the immediate ancestor
	///
	/// @return Pointer to immediate ancestor or null if there's no ancestor.
    ConstGenParticlePtr ancestor(ConstGenParticlePtr p) const;

    /// Get the immediate ancestor of a particle, which is @e not an
	/// ignored particle.
	///
	/// @param p Particle for which to get the immediate ancestor
	///
	/// @return Pointer to immediate ancestor or null if there's no ancestor.
    ConstGenParticlePtr ancestor(ConstGenParticlePtr p, bool) const;

    /// Particle types to test for
    std::vector<int> _pdgIds;

  };


}

#endif
