// -*- C++ -*-
#ifndef RIVET_ParticleFinder_HH
#define RIVET_ParticleFinder_HH

#include "Rivet/Projection.hh"
#include "Rivet/Cuts.hh"

namespace Rivet {


  /// @brief Base class for projections which return subsets of an event's particles
  class ParticleFinder : public Projection {
  public:

    /// @name Object lifetime management
    //@{

    /// Construction using Cuts object
    ParticleFinder(const Cut& c=Cuts::open())
      : _cuts(c), _theParticles()
    { }

    /// Virtual destructor for inheritance
    virtual ~ParticleFinder() {}

    /// Clone on the heap.
    virtual const Projection* clone() const = 0;

    //@}


    /// @name Particle accessors
    //@{

    /// Get the final-state particles in no particular order, with no cuts.
    virtual const Particles& particles() const { return _theParticles; }

    /// Access the projected final-state particles.
    size_t size() const { return particles().size(); }

    /// Is this final state empty?
    bool empty() const { return particles().empty(); }
    /// @deprecated Is this final state empty?
    DEPRECATED("Use empty()")
    bool isEmpty() const { return particles().empty(); }


    /// @brief Get the final-state particles, with optional cuts.
    /// @note Returns a copy rather than a reference, due to cuts
    /// @todo Can't this be a const Cut& arg?
    Particles particles(const Cut & c) const {
      // Just return a copy of particles() if the cut is open
      if (c == Cuts::open()) return particles();
      // If there is a non-trivial cut...
      Particles rtn;
      rtn.reserve(size());
      foreach (const Particle& p, particles())
        if (c->accept(p)) rtn.push_back(p);
      return rtn;
    }

    /// Get the final-state particles, ordered by supplied sorting function object.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @todo Can't this be a const Cut& arg?
    /// @todo Use a std::function instead of typename F?
    template <typename F>
    Particles particles(F sorter, const Cut & c=Cuts::open()) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(particles(c), sorter);
    }

    /// Get the final-state particles, ordered by supplied sorting function object.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @todo Can't this be a const Cut& arg?
    /// @todo Use a std::function instead of typename F?
    template <typename F>
    Particles particles(const Cut & c, F sorter) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(particles(c), sorter);
    }

    /// Get the final-state particles, ordered by decreasing \f$ p_T \f$ and with optional cuts.
    ///
    /// This is a very common use-case, so is available as syntatic sugar for particles(c, cmpMomByPt).
    Particles particlesByPt(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByPt);
    }

    /// Get the final-state particles, ordered by decreasing \f$ p_T \f$ and with a cut on minimum \f$ p_T \f$.
    ///
    /// This is a very common use-case, so is available as syntatic sugar for particles(Cuts::pT >= ptmin, cmpMomByPt).
    Particles particlesByPt(double ptmin) const {
      return particles(Cuts::pT >= ptmin, cmpMomByPt);
    }


    /// @name Little-used sorted accessors
    /// @deprecated Use the versions with a sorter function argument
    //@{

    /// Get the final-state particles, ordered by decreasing \f$ p \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByP(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByP);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByE(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByE);
    }

    /// Get the final-state particles, ordered by decreasing \f$ E_T \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByEt(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByEt);
    }

    /// Get the final-state particles, ordered by increasing \f$ \eta \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByEta(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByEta);
    }

    /// Get the final-state particles, ordered by increasing \f$ |\eta| \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByModEta(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByAbsEta);
    }

    /// Get the final-state particles, ordered by increasing \f$ y \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByRapidity(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByRap);
    }

    /// Get the final-state particles, ordered by increasing \f$ |y| \f$.
    /// @todo Remove, since there is the templated method or sortByX methods available for these unusual cases?
    /// @deprecated Use the version with a sorter function argument
    DEPRECATED("Use the version with a sorter function argument")
    Particles particlesByModRapidity(const Cut & c=Cuts::open()) const {
      return particles(c, cmpMomByAbsRap);
    }

    //@}

    //@}


    /// @todo Replace with cuts() accessor
    ///virtual Cut cuts() const { return _cuts; }


    /// @name For JetAlg compatibility
    //@{

    typedef Particle entity_type;
    typedef Particles collection_type;

    /// Template-usable interface common to JetAlg
    const collection_type& entities() const {
      return particles();
    }

    //@}


  protected:

    /// Apply the projection to the event
    virtual void project(const Event& e) = 0;

    /// Compare projections
    virtual int compare(const Projection& p) const;

    /// The kinematic cuts cuts
    Cut _cuts;

    /// The found particles returned by the particles() methods
    Particles _theParticles;

  };


}

#endif
