// -*- C++ -*-
#ifndef RIVET_Projection_HH
#define RIVET_Projection_HH

#include "Rivet/Projection.fhh"
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Cuts.hh"
// NOTE: Cmp.hh, Event.hh and Particle.hh included at the bottom

namespace Rivet {


  // Forward declaration
  class Event;


  /// @brief Base class for all Rivet projections.
  ///
  /// Projection is the base class of all Projections to be used by
  /// Rivet. A Projection object can be assigned to an Event object and
  /// will then define a processed part of the information available in
  /// the Event, which then can be used by other Projection objects
  /// and/or Analysis objects.
  ///
  /// The main virtual functions to be overridden by concrete sub-classes
  /// are project(const Event &) and compare(const Projection &).
  class Projection : public ProjectionApplier {
  public:

    /// Event is a friend.
    friend class Event;

    /// The Cmp specialization for Projection is a friend.
    friend class Cmp<Projection>;



    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor.
    Projection();

    /// Clone on the heap.
    virtual unique_ptr<Projection> clone() const = 0;

    /// The destructor.
    virtual ~Projection();

    //@}


    /// Get the name of the projection.
    virtual std::string name() const {
      return _name;
    }

    /// Get the state of the projetion.
    bool valid() const {
      return _isValid;
    }

    /// Get the state of the projetion.
    bool failed() const {
      return !valid();
    }

    /// @name Projection operation and comparison
    //@{

    /// Take the information available in the Event and make the
    /// calculations necessary to obtain the projection. Note that this
    /// function must never be called except inside the
    /// Event::applyProjection(Projection *) function.
    virtual void project(const Event& e) = 0;

    /// This function is used to define a unique ordering between
    /// different Projection objects of the same class. If this is
    /// considered to be equivalent to the Projector object, \a p, in the
    /// argument the function should return 0. If this object should be
    /// ordered before \a p a negative value should be returned,
    /// otherwise a positive value should be returned. This function must
    /// never be called explicitly, but should only be called from the
    /// operator<(const Projection &). When implementing the function in
    /// concrete sub-classes, it is then guaranteed that the Projection
    /// object \a p in the argument is of the same class as the sub-class
    /// and can be safely dynamically casted to that class.
    ///
    /// When implementing this function in a sub-class, the immediate
    /// base class version of the function should be called first. If the
    /// base class function returns a non-zero value, that value should
    /// be returned immediately. Only if zero is returned should this
    /// function check the member variables of the sub-class to determine
    /// whether this should be ordered before or after \a p, or if it is
    /// equivalent with \a p.
    virtual int compare(const Projection& p) const = 0;

    /// Determine whether this object should be ordered before the object
    /// \a p given as argument. If \a p is of a different class than
    /// this, the before() function of the corresponding type_info
    /// objects is used. Otherwise, if the objects are of the same class,
    /// the virtual compare(const Projection &) will be returned.
    bool before(const Projection& p) const;

    //@}



    /// @name Beam configuration
    /// @todo Does it really make sense to restrict Projections to particular beam configs? Do we use this in practice?
    //@{

    /// Return the allowed beam pairs on which this projection can operate, not
    /// including recursion. Derived classes should ensure that all contained
    /// projections are registered in the @a _projections set for the beam
    /// constraint chaining to work.
    /// @todo Remove the beam constraints system from projections.
    virtual const std::set<PdgIdPair> beamPairs() const;


    /// Add a colliding beam pair.
    /// @todo This deserves a better name!
    Projection& addPdgIdPair(PdgId beam1, PdgId beam2) {
      _beamPairs.insert(PdgIdPair(beam1, beam2));
      return *this;
    }

    //@}


  protected:

    /// Get a Log object based on the getName() property of the calling projection object.
    Log& getLog() const {
      string logname = "Rivet.Projection." + name();
      return Log::getLog(logname);
    }

    /// Used by derived classes to set their name.
    void setName(const std::string& name) {
      _name = name;
    }

    /// Set the projection in an unvalid state.
    void fail() {
      _isValid = false;
    }

    /// Shortcut to make a named Cmp<Projection> comparison with the @c *this
    /// object automatically passed as one of the parent projections.
    Cmp<Projection> mkNamedPCmp(const Projection& otherparent, const std::string& pname) const;

    /// Shortcut to make a named Cmp<Projection> comparison with the @c *this
    /// object automatically passed as one of the parent projections.
    ///
    /// @note Alias for mkNamedPCmp
    Cmp<Projection> mkPCmp(const Projection& otherparent, const std::string& pname) const;

    /// Block Projection copying
    virtual Projection& operator = (const Projection&);


  private:

    /// Name variable is used by the base class messages to identify
    /// which derived class is being handled.
    string _name;

    /// Beam-type constraint.
    /// @todo Remove?
    set<PdgIdPair> _beamPairs;

    /// Flag to tell if this projection is in a valid state.
    bool _isValid;
    
  };


}


/// Define "less" operator for Projection* containers in terms of the Projection::before virtual method.
inline bool std::less<const Rivet::Projection *>::operator()(const Rivet::Projection* x,
                                                             const Rivet::Projection* y) const {
  return x->before(*y);
}


#endif


#include "Rivet/Event.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Tools/Cmp.hh"


/// @def DEFAULT_RIVET_PROJ_CLONE
/// Preprocessor define to prettify the manky constructor with name string argument
#define DEFAULT_RIVET_PROJ_CLONE(clsname) \
  virtual unique_ptr<Projection> clone() const { return unique_ptr<Projection>(new clsname(*this)); }
