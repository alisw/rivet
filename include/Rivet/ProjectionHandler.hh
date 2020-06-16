// -*- C++ -*-
#ifndef RIVET_ProjectionHandler_HH
#define RIVET_ProjectionHandler_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Projection.fhh"

namespace Rivet {


  /// Typedef for Projection (smart) pointer
  typedef std::shared_ptr<const Projection> ProjHandle;

  // Forward declaration.
  class ProjectionApplier;

  /// @brief The projection handler is a central repository for projections to be used
  /// in a Rivet analysis run.
  ///
  /// Without centralised projections, it can be hard to know which of an
  /// equivalent set of projections will be run on a particular event. In turn,
  /// this may mean that certain projections in the chain can go out of scope
  /// unexpectedly. There were originally also the issues that projections may
  /// need to be held as member pointers to an abstract base class, since
  /// post-construction setup is needed; that projections contained pointers to
  /// their own dependency chain, which could go out of scope; and that
  /// projection members could be modified after being applied to an event
  /// which, due to the caching model, would have unpredictable consequences.
  ///
  /// By centralising all the projections, these issues are eliminated, as well
  /// as allowing analysis classes to contain fewer data members (since
  /// projections are now better accessed by name than by storing a data member
  /// reference or pointer).
  ///
  /// The core of the ProjectionHandler design is that it is a singleton class,
  /// essentially a wrapper around a map of @c Projection*, indexed by a hash of
  /// the registering object and its local name for the registered projection.
  ///
  class ProjectionHandler {
  public:

    /// ProjectionApplier's destructor needs to trigger cleaning up the proj handler repo
    friend class ProjectionApplier;

    /// Typedef for a vector of Projection pointers.
    typedef set<ProjHandle> ProjHandles;

    /// @brief Typedef for the structure used to contain named projections for a
    /// particular containing Analysis or Projection.
    typedef map<string, ProjHandle> NamedProjs;

    /// Enum to specify depth of projection search.
    enum ProjDepth { SHALLOW, DEEP };


  private:

    /// Structure used to map a containing Analysis or Projection to its set of
    /// contained projections.
    typedef map<const ProjectionApplier*, NamedProjs> NamedProjsMap;

    /// Core data member, associating a given containing class (via a
    /// ProjectionApplier pointer) to its contained projections.
    NamedProjsMap _namedprojs;

    /// Cache of {@link Projection}s for reverse lookup, to speed up registering
    /// new projections as @c _namedprojs gets large.
    ProjHandles _projs;


  private:

    /// @name Construction. */
    //@{

    /// Private destructor means no inheritance from this class.
    ~ProjectionHandler() = default;

    /// The assignment operator is hidden.
    ProjectionHandler& operator = (const ProjectionHandler&) = delete;

    /// The copy constructor is hidden.
    ProjectionHandler(const ProjectionHandler&) = delete;

    /// The standard constructor.
    ProjectionHandler() = default;

    //@}


  public:

    /// Singleton creation function
    static ProjectionHandler& getInstance() {
      static ProjectionHandler _instance;
      return _instance;
    }


  public:

    /// @name Projection registration
    //@{
    /// Attach and retrieve a projection as a reference.
    const Projection& registerProjection(const ProjectionApplier& parent,
                                         const Projection& proj,
                                         const string& name);
    //@}


  private:

    /// @name Projection registration internal helpers
    //@{

    /// Try to get an equivalent projection from the system
    /// @returns 0 if no equivalent projection found
    ProjHandle _getEquiv(const Projection& proj) const;

    /// Make a clone of proj, copying across child references from the original
    unique_ptr<Projection> _clone(const Projection& proj);

    /// Internal function to do the registering
    const Projection& _register(const ProjectionApplier& parent,
                                ProjHandle proj,
                                const string& name);

    /// Get a string dump of the current ProjHandler structure
    string _getStatus() const;

    /// Check that this parent projection doesn't already use this name
    bool _checkDuplicate(const ProjectionApplier& parent,
                         const Projection& proj,
                         const string& name) const;

    //@}


  public:

    /// @name Projection retrieval. */
    //@{

    /// Check if there is a @a name projection registered by @a parent
    bool hasProjection(const ProjectionApplier& parent, const string& name) const;

    /// Retrieve a named projection for the given parent. Returning as a
    /// reference is partly to discourage ProjectionApplier classes from storing
    /// pointer members to the registered projections, since that can lead to
    /// problems and there is no need to do so.
    const Projection& getProjection(const ProjectionApplier& parent,
                                    const string& name) const;

    /// Get child projections for the given parent. By default this will just
    /// return the projections directly contained by the @a parent, but the @a
    /// depth argument can be changed to do a deep retrieval, which will recurse
    /// through the whole projection chain. In this case, there is no protection
    /// against getting stuck in a circular projection dependency loop.
    set<const Projection*> getChildProjections(const ProjectionApplier& parent,
                                               ProjDepth depth=SHALLOW) const;
    //@}


  private:

    /// Remove a ProjectionApplier: designed to only be called by ~ProjectionApplier (as a friend)
    void removeProjectionApplier(ProjectionApplier& parent);


  private:

    // /// Get map of named projections belonging to @a parent.
    // /// Throws an exception if @a parent has not got any registered projections.
    // const NamedProjs& namedProjs(const ProjectionApplier* parent) const {
    //   NamedProjsMap::const_iterator nps = _namedprojs.find(parent);
    //   if (nps == _namedprojs.end()) {
    //     stringstream ss;
    //     ss << "No NamedProjs registered for parent " << parent;
    //     throw Error(ss.str());
    //   }
    //   return *nps;
    // }


  };


}

#endif
