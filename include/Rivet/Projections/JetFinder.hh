// -*- C++ -*-
#ifndef RIVET_JetFinder_HH
#define RIVET_JetFinder_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  /// Abstract base class for projections which can return a set of {@link Jet}s.
  class JetFinder : public Projection {
  public:

    /// Enum for the treatment of muons: whether to include all, some, or none in jet-finding
    enum class Muons { NONE, DECAY, ALL };

    /// Enum for the treatment of invisible particles: whether to include all, some, or none in jet-finding
    enum class Invisibles { NONE, DECAY, ALL };



    /// Constructor
    JetFinder(const FinalState& fs,
           Muons usemuons = Muons::ALL,
           Invisibles useinvis = Invisibles::NONE);

    /// Default constructor
    JetFinder() = default;

    /// Clone on the heap.
    virtual unique_ptr<Projection> clone() const = 0;

    /// Destructor
    virtual ~JetFinder() = default;


    /// @name Control the treatment of muons and invisible particles
    ///
    /// Since MC-based jet calibration (and/or particle flow) can add back in
    /// particles that weren't seen in calorimeters/trackers.
    //@{

    /// @brief Include (some) muons in jet construction.
    ///
    /// The default behaviour is that jets are only constructed from visible
    /// particles. Some jet studies, including those from ATLAS, use a definition
    /// in which neutrinos from hadron decays are included via MC-based calibrations.
    /// Setting this flag to true avoids the automatic restriction to a VisibleFinalState.
    void useMuons(Muons usemuons = Muons::ALL) {
      _useMuons = usemuons;
    }

    /// @brief Include (some) invisible particles in jet construction.
    ///
    /// The default behaviour is that jets are only constructed from visible
    /// particles. Some jet studies, including those from ATLAS, use a definition
    /// in which neutrinos from hadron decays are included via MC-based calibrations.
    /// Setting this flag to true avoids the automatic restriction to a VisibleFinalState.
    void useInvisibles(Invisibles useinvis = Invisibles::DECAY) {
      _useInvisibles = useinvis;
    }

    /// @brief obsolete chooser
    DEPRECATED("make an explicit choice from Invisibles::{NONE,DECAY,ALL}. This boolean call does not allow for ALL")
    void useInvisibles(bool useinvis) {
      _useInvisibles = useinvis ? Invisibles::DECAY : Invisibles::NONE;
    }

    //@}


    /// @name Access to jet objects
    //@{

    /// Get jets in no guaranteed order, with an optional Cut
    /// @note Returns a copy rather than a reference, due to cuts
    virtual Jets jets(const Cut& c=Cuts::open()) const {
      return filter_select(_jets(), c);
    }

    /// Get jets in no guaranteed order, with a selection functor
    /// @note Returns a copy rather than a reference, due to cuts
    virtual Jets jets(const JetSelector& selector) const {
      return filter_select(_jets(), selector);
    }


    /// Get the jets with a Cut applied, and ordered by supplied sorting functor
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    Jets jets(const Cut& c, const JetSorter& sorter) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(jets(c), sorter);
    }

    /// Get the jets, ordered by supplied sorting functor, with an optional Cut
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    Jets jets(const JetSorter& sorter, const Cut& c=Cuts::open()) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return jets(c, sorter);
    }

    /// Get the jets, ordered by supplied sorting function object, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    Jets jets(const JetSelector& selector, const JetSorter& sorter) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(jets(selector), sorter);
    }

    /// Get the jets, ordered by supplied sorting functor and with a selection functor applied
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    Jets jets(const JetSorter& sorter, const JetSelector selector) const {
      return jets(selector, sorter);
    }


    /// Get the jets, ordered by \f$ p_T \f$, with optional cuts.
    ///
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(c, cmpMomByPt).
    Jets jetsByPt(const Cut& c=Cuts::open()) const {
      return jets(c, cmpMomByPt);
    }

    /// Get the jets, ordered by \f$ p_T \f$, with cuts via a selection functor.
    ///
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(c, cmpMomByPt).
    Jets jetsByPt(const JetSelector& selector) const {
      return jets(selector, cmpMomByPt);
    }

    /// Get the jets, ordered by \f$ p_T \f$, with a cut on \f$ p_\perp \f$.
    ///
    /// @deprecated Use the version with a Cut argument
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(Cuts::pT >= ptmin, cmpMomByPt).
    Jets jetsByPt(double ptmin) const {
      return jets(Cuts::pT >= ptmin, cmpMomByPt);
    }

    //@}


  protected:

    /// @brief Internal pure virtual method for getting jets in no guaranteed order.
    virtual Jets _jets() const = 0;


  public:

    /// Count the jets
    size_t size() const { return jets().size(); }
    /// Count the jets after a Cut is applied.
    size_t size(const Cut& c) const { return jets(c).size(); }
    /// Count the jets after a selection functor is applied.
    size_t size(const JetSelector& s) const { return jets(s).size(); }

    /// Is this jet finder empty?
    bool empty() const { return size() == 0; }
    /// Is this jet finder empty after a Cut is applied?
    bool empty(const Cut& c) const { return size(c) == 0; }
    /// Is this jet finder empty after a selection functor is applied?
    bool empty(const JetSelector& s) const { return size(s) == 0; }

    /// Clear the projection.
    virtual void reset() = 0;

    typedef Jet entity_type;
    typedef Jets collection_type;

    /// Template-usable interface common to FinalState.
    collection_type entities() const { return jets(); }

    // /// Do the calculation locally (no caching).
    // virtual void calc(const Particles& constituents, const Particles& tagparticles=Particles()) = 0;


  protected:

    /// Perform the projection on the Event.
    virtual void project(const Event& e) = 0;

    /// Compare projections.
    virtual CmpState compare(const Projection& p) const = 0;


  protected:

    /// Flag to determine whether or not to exclude (some) muons from the would-be constituents.
    Muons _useMuons;

    /// Flag to determine whether or not to exclude (some) invisible particles from the would-be constituents.
    Invisibles _useInvisibles;


  };


}

#endif
