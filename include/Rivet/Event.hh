// -*- C++ -*-
#ifndef RIVET_Event_HH
#define RIVET_Event_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"

namespace Rivet {


  /// Rivet wrapper for HepMC event and Projection references.
  ///
  /// Event is a concrete class representing an generated event in Rivet. It is
  /// constructed given a HepMC::GenEvent, a pointer to which is kept by the
  /// Event object throughout its lifetime. The user must therefore make sure
  /// that the corresponding HepMC::GenEvent will persist at least as long as
  /// the Event object.
  ///
  /// In addition to the HepMC::GenEvent object the Event also keeps track of
  /// all Projection objects which have been applied to the Event so far.
  class Event {
  public:

    /// @name Constructors and destructors.
    //@{

    /// Constructor from a HepMC GenEvent pointer
    Event(const GenEvent* ge)
      : _genevent_original(ge), _genevent(*ge)
    { assert(ge); _init(*ge); }

    /// Constructor from a HepMC GenEvent reference
    /// @deprecated HepMC uses pointers, so we should talk to HepMC via pointers
    Event(const GenEvent& ge)
      : _genevent_original(&ge), _genevent(ge)
    { _init(ge); }

    /// Copy constructor
    Event(const Event& e)
      : _genevent_original(e._genevent_original), _genevent(e._genevent)
    {  }

    //@}


    /// @name Major event properties
    //@{

    /// The generated event obtained from an external event generator
    const GenEvent* genEvent() const { return &_genevent; }

    /// @brief The generation weight associated with the event
    ///
    /// @todo This needs to be revisited when we finally add the mechanism to
    /// support NLO counter-events and weight vectors.
    double weight() const;

    /// Get the beam particles
    ParticlePair beams() const;

    /// Get the beam centre-of-mass energy
    double sqrtS() const;

    /// Get the beam centre-of-mass energy per nucleon
    double asqrtS() const;

    /// Get the generator centrality (impact-parameter quantile in [0,1]; or -1 if undefined (usual for non-HI generators))
    double centrality() const;

    // /// Get the boost to the beam centre-of-mass
    // Vector3 beamCMSBoost() const;

    // /// Get the boost to the beam centre-of-mass
    // LorentzTransform beamCMSTransform();

    //@}


    /// @name Access to event particles
    //@{

    /// All the raw GenEvent particles, wrapped in Rivet::Particle objects
    const Particles& allParticles() const;

    /// @brief All the raw GenEvent particles, wrapped in Rivet::Particle objects, but with a Cut applied
    ///
    /// @note Due to the cut, this returns by value, i.e. involves an expensive copy
    inline Particles allParticles(const Cut& c) const {
      return filter_select(allParticles(), c);
    }

    /// @brief All the raw GenEvent particles, wrapped in Rivet::Particle objects, but with a selection function applied
    ///
    /// @note Due to the cut, this returns by value, i.e. involves an expensive copy
    template <typename FN>
    inline Particles allParticles(const FN& f) const {
      return filter_select(allParticles(), f);
    }

    //@}


    /// @name Projection running
    //@{

    /// @brief Add a projection @a p to this Event.
    ///
    /// If an equivalent Projection has been applied before, the
    /// Projection::project(const Event&) of @a p is not called and a reference
    /// to the previous equivalent projection is returned. If no previous
    /// Projection was found, the Projection::project(const Event&) of @a p is
    /// called and a reference to @a p is returned.
    template <typename PROJ>
    const PROJ& applyProjection(PROJ& p) const {
      Log& log = Log::getLog("Rivet.Event");
      log << Log::TRACE << "Applying projection " << &p << " (" << p.name() << ") -> comparing to projections " << _projections << endl;
      // First search for this projection *or an equivalent* in the already-executed list
      const Projection* cpp(&p);
      std::set<const Projection*>::const_iterator old = _projections.find(cpp);
      if (old != _projections.end()) {
        log << Log::TRACE << "Equivalent projection found -> returning already-run projection " << *old << endl;
        const Projection& pRef = **old;
        return pcast<PROJ>(pRef);
      }
      // If this one hasn't been run yet on this event, run it and add to the list
      log << Log::TRACE << "No equivalent projection in the already-run list -> projecting now" << endl;
      Projection* pp = const_cast<Projection*>(cpp);
      pp->project(*this);
      _projections.insert(pp);
      return p;
    }


    /// @brief Add a projection @a p to this Event by pointer.
    template <typename PROJ>
    const PROJ& applyProjection(PROJ* pp) const {
      if (!pp) throw Error("Event::applyProjection(PROJ*): Projection pointer is null.");
      return applyProjection(*pp);
    }

    //@}


  private:

    /// @brief Actual (shared) implementation of the constructors from GenEvents
    void _init(const GenEvent& ge);

    // /// @brief Convert the GenEvent to use conventional alignment
    // ///
    // /// For example, FHerwig only produces DIS events in the unconventional
    // /// hadron-lepton orientation and has to be corrected for DIS analysis
    // /// portability.
    // void _geNormAlignment();

    /// @brief The generated event, as obtained from an external generator.
    ///
    /// This is the original GenEvent. In practise the version seen by users
    /// will often/always be a modified one.
    ///
    /// @todo Provide access to this via an Event::originalGenEvent() method? If requested...
    const GenEvent* _genevent_original;

    /// @brief The GenEvent used by Rivet analysis projections etc.
    ///
    /// This version may be rotated to a "normal" alignment, have
    /// generator-specific particles stripped out, etc.  If an analysis is
    /// affected by these modifications, it is probably an unphysical analysis!
    ///
    /// Stored as a non-pointer since it may get overwritten, and memory for
    /// copying and cleanup is neater this way.
    /// @todo Change needed for HepMC3?
    mutable GenEvent _genevent;

    /// All the GenEvent particles, wrapped as Rivet::Particles
    /// @note To be populated lazily, hence mutability
    mutable Particles _particles;

    /// The set of Projection objects applied so far
    mutable std::set<ConstProjectionPtr> _projections;

  };


}

#endif
