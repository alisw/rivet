// -*- C++ -*-
#ifndef RIVET_JetAlg_HH
#define RIVET_JetAlg_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  /// Abstract base class for projections which can return a set of {@link Jet}s.
  class JetAlg : public Projection {
  public:

    /// Enum for the treatment of muons: whether to include all, some, or none in jet-finding
    enum MuonsStrategy { NO_MUONS, DECAY_MUONS, ALL_MUONS };

    /// Enum for the treatment of invisible particles: whether to include all, some, or none in jet-finding
    enum InvisiblesStrategy { NO_INVISIBLES, DECAY_INVISIBLES, ALL_INVISIBLES };



    /// Constructor
    JetAlg(const FinalState& fs, MuonsStrategy usemuons=JetAlg::ALL_MUONS, InvisiblesStrategy useinvis=JetAlg::NO_INVISIBLES);

    /// Default constructor
    JetAlg() {};

    /// Clone on the heap.
    virtual unique_ptr<Projection> clone() const = 0;

    /// Destructor
    virtual ~JetAlg() { }


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
    void useMuons(MuonsStrategy usemuons=ALL_MUONS) {
      _useMuons = usemuons;
    }

    /// @brief Include (some) invisible particles in jet construction.
    ///
    /// The default behaviour is that jets are only constructed from visible
    /// particles. Some jet studies, including those from ATLAS, use a definition
    /// in which neutrinos from hadron decays are included via MC-based calibrations.
    /// Setting this flag to true avoids the automatic restriction to a VisibleFinalState.
    void useInvisibles(InvisiblesStrategy useinvis=DECAY_INVISIBLES) {
      _useInvisibles = useinvis;
    }

    /// @brief Include (some) invisible particles in jet construction.
    ///
    /// The default behaviour is that jets are only constructed from visible
    /// particles. Some jet studies, including those from ATLAS, use a definition
    /// in which neutrinos from hadron decays are included via MC-based calibrations.
    /// Setting this flag to true avoids the automatic restriction to a VisibleFinalState.
    ///
    /// @deprecated Use the enum-arg version instead. Will be removed in Rivet v3
    void useInvisibles(bool useinvis) {
      _useInvisibles = useinvis ? DECAY_INVISIBLES : NO_INVISIBLES;
    }

    //@}


    /// @name Access to jet objects
    //@{

    /// Get jets in no guaranteed order, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts
    virtual Jets jets(const Cut& c=Cuts::open()) const {
      return filterBy(_jets(), c);
      // const Jets rawjets = _jets();
      // // Just return a copy of rawjets if the cut is open
      // if (c == Cuts::open()) return rawjets;
      // // If there is a non-trivial cut...
      // /// @todo Use an STL erase(remove_if) and lambda function for this
      // Jets rtn;
      // rtn.reserve(size());
      // foreach (const Jet& j, rawjets)
      //   if (c->accept(j)) rtn.push_back(j);
      // return rtn;
    }

    /// Get the jets, ordered by supplied sorting function object, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    template <typename F>
    Jets jets(F sorter, const Cut& c=Cuts::open()) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(jets(c), sorter);
    }

    /// Get the jets, ordered by supplied sorting function object, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    template <typename F>
    Jets jets(const Cut& c, F sorter) const {
      /// @todo Will the vector be efficiently std::move'd by value through this function chain?
      return sortBy(jets(c), sorter);
    }


    /// Get the jets, ordered by \f$ p_T \f$, with optional cuts.
    ///
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(c, cmpMomByPt).
    /// @todo The other sorted accessors should be removed in a cleanup.
    Jets jetsByPt(const Cut& c=Cuts::open()) const {
      return jets(c, cmpMomByPt);
    }

    //@}


    /// @name Old sorted jet accessors
    /// @deprecated Use the versions with sorter function arguments. These will be removed in Rivet v3
    //@{

    /// Get the jets, ordered by \f$ |p| \f$, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @deprecated Use the version with a sorter function argument.
    DEPRECATED("Use the version with a sorter function argument.")
    Jets jetsByP(const Cut& c=Cuts::open()) const {
      return jets(c, cmpMomByP);
    }

    /// Get the jets, ordered by \f$ E \f$, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @deprecated Use the version with a sorter function argument.
    DEPRECATED("Use the version with a sorter function argument.")
    Jets jetsByE(const Cut &c=Cuts::open()) const {
      return jets(c, cmpMomByE);
    }

    /// Get the jets, ordered by \f$ E_T \f$, with optional cuts on \f$ p_\perp \f$ and rapidity.
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    /// @deprecated Use the version with a sorter function argument.
    DEPRECATED("Use the version with a sorter function argument.")
    Jets jetsByEt(const Cut& c=Cuts::open()) const {
      return jets(c, cmpMomByEt);
    }

    //@}


    /// @name Old jet accessors
    /// @deprecated Use the versions with Cut arguments
    //@{

    /// Get jets in no guaranteed order, with optional cuts on \f$ p_\perp \f$ and rapidity.
    ///
    /// @deprecated Use the version with a Cut argument
    /// @note Returns a copy rather than a reference, due to cuts
    DEPRECATED("Use the version with a Cut argument.")
    Jets jets(double ptmin, double ptmax=MAXDOUBLE,
              double rapmin=-MAXDOUBLE, double rapmax=MAXDOUBLE,
              RapScheme rapscheme=PSEUDORAPIDITY) const {
      if (rapscheme == PSEUDORAPIDITY) {
        return jets((Cuts::pT >= ptmin) & (Cuts::pT < ptmax) & (Cuts::rapIn(rapmin, rapmax)));
      } else if (rapscheme == RAPIDITY) {
        return jets((Cuts::pT >= ptmin) & (Cuts::pT < ptmax) & (Cuts::etaIn(rapmin, rapmax)));
      }
      throw LogicError("Unknown rapidity scheme. This shouldn't be possible!");
    }

    /// Get the jets, ordered by \f$ p_T \f$, with a cut on \f$ p_\perp \f$.
    ///
    /// @deprecated Use the version with a Cut argument
    /// @note Returns a copy rather than a reference, due to cuts and sorting
    ///
    /// This is a very common use-case, so is available as syntatic sugar for jets(Cuts::pT >= ptmin, cmpMomByPt).
    /// @todo The other sorted accessors should be removed in a cleanup.
    Jets jetsByPt(double ptmin) const {
      return jets(Cuts::pT >= ptmin, cmpMomByPt);
    }

    //@}


  protected:

    /// @brief Internal pure virtual method for getting jets in no guaranteed order.
    virtual Jets _jets() const = 0;


  public:

    /// Number of jets passing the provided Cut.
    size_t numJets(const Cut& c=Cuts::open()) const { return jets(c).size(); }

    /// Number of jets (without cuts).
    size_t size() const { return jets().size(); }
    /// Whether the inclusive jet collection is empty.
    bool empty() const { return size() != 0; }

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
    virtual int compare(const Projection& p) const = 0;


  protected:

    /// Flag to determine whether or not to exclude (some) muons from the would-be constituents.
    MuonsStrategy _useMuons;

    /// Flag to determine whether or not to exclude (some) invisible particles from the would-be constituents.
    InvisiblesStrategy _useInvisibles;


  };


}

#endif
