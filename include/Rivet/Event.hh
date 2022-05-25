// -*- C++ -*-
#ifndef RIVET_Event_HH
#define RIVET_Event_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"

namespace Rivet {


  /// @brief Representation of a HepMC event, and enabler of Projection caching
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
    Event(const GenEvent* ge, const vector<size_t>& indices = {}, bool strip = false)
      : _weightIndices(indices), _genevent_original(ge) {
      assert(ge);
      _genevent = *ge;
      if ( strip ) _strip(_genevent);
      _init(*ge);
    }

    /// Constructor from a HepMC GenEvent reference
    /// @deprecated HepMC uses pointers, so we should talk to HepMC via pointers
    Event(const GenEvent& ge, const vector<size_t>& indices = {}, bool strip = false)
      : _weightIndices(indices), _genevent_original(&ge), _genevent(ge) {
        if ( strip ) _strip(_genevent);
        _init(ge);
      }

    /// Copy constructor
    Event(const Event& e)
      : _weightIndices(e._weightIndices),
        _genevent_original(e._genevent_original),
        _genevent(e._genevent)
    {  }

    //@}


    /// @name Major event properties
    //@{

    /// The generated event obtained from an external event generator
    const GenEvent* genEvent() const { return &_genevent; }

    /// The generated event obtained from an external event generator
    const GenEvent* originalGenEvent() const { return _genevent_original; }

    /// Get the beam particles
    ParticlePair beams() const;

    /// Get the beam centre-of-mass energy
    double sqrtS() const;

    /// Get the beam centre-of-mass energy per nucleon
    double asqrtS() const;

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

    /// @brief The generation weights associated with the event
    std::valarray<double> weights() const;

    /// @brief The generation cross-sections associated with the event
    std::vector<std::pair<double, double>> crossSections() const;

    /// @brief Obsolete weight method. Always returns 1 now.
    DEPRECATED("Event weight does not need to be included anymore. For compatibility, it's always == 1 now.")
    double weight() const { return 1.0; }
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
    ///
    /// @todo Can make this non-templated, since only cares about ptr to Projection base class
    ///
    /// @note Comparisons here are by direct pointer comparison, because
    /// equivalence is guaranteed if pointers are equal, and inequivalence
    /// guaranteed if they aren't, thanks to the ProjectionHandler registry
    template <typename PROJ>
    const PROJ& applyProjection(PROJ& p) const {
      static bool docaching = getEnvParam("RIVET_CACHE_PROJECTIONS", true);
      if (docaching) {
        MSG_TRACE("Applying projection " << &p << " (" << p.name() << ") -> comparing to projections " << _projections);
        // First search for this projection *or an equivalent* in the already-executed list
        const Projection* cpp(&p);
        /// @note Currently using reint cast to integer type to bypass operator==(Proj*, Proj*)
        // std::set<const Projection*>::const_iterator old = _projections.find(cpp);
        std::set<const Projection*>::const_iterator old = std::begin(_projections);
        std::uintptr_t recpp = reinterpret_cast<std::uintptr_t>(cpp);
        for (; old != _projections.end(); ++old)
          if (reinterpret_cast<std::uintptr_t>(*old) == recpp) break;
        if (old != _projections.end()) {
          MSG_TRACE("Equivalent projection found -> returning already-run projection " << *old);
          const Projection& pRef = **old;
          return pcast<PROJ>(pRef);
        }
        MSG_TRACE("No equivalent projection in the already-run list -> projecting now");
      } else {
        MSG_TRACE("Applying projection " << &p << " (" << p.name() << ") WITHOUT projection caching & comparison");
      }
      // If this one hasn't been run yet on this event, run it and add to the list
      Projection* pp = const_cast<Projection*>(&p);
      pp->_isValid = true;
      pp->project(*this);
      if (docaching) _projections.insert(pp);
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

    /// Get a Log object for Event
    Log& getLog() const;

    /// @brief Actual (shared) implementation of the constructors from GenEvents
    void _init(const GenEvent& ge);

    /// @brief Remove uninteresting or unphysical particles in the
    /// GenEvent to speed up searches.
    void _strip(GenEvent & ge);

    // /// @brief Convert the GenEvent to use conventional alignment
    // ///
    // /// For example, FHerwig only produces DIS events in the unconventional
    // /// hadron-lepton orientation and has to be corrected for DIS analysis
    // /// portability.
    // void _geNormAlignment();

    /// @brief The indices of the selected weights, as instructed to the AnalysisHandler.
    ///
    /// The user can (de-)select weights and the AnalysisHandler knows about the subset
    /// that match the specifications.
    const std::vector<size_t> _weightIndices;

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

    /// A cached set of event weights (a single unit weight if the original weight vector was empty)
    mutable std::valarray<double> _weights;

    /// A cached set of event cross-section & errors (a pair of zeros if the original weight vector was empty)
    mutable std::vector<std::pair<double,double>> _xsecs;
  };


}

#endif
