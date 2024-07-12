#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

#include "Rivet/Config/RivetCommon.hh"
#include "YODA/AnalysisObject.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

#include <map>
#include <valarray>

namespace YODA {

  typedef std::shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;

  typedef std::shared_ptr<YODA::Counter> CounterPtr;
  typedef std::shared_ptr<YODA::Histo1D> Histo1DPtr;
  typedef std::shared_ptr<YODA::Histo2D> Histo2DPtr;
  typedef std::shared_ptr<YODA::Profile1D> Profile1DPtr;
  typedef std::shared_ptr<YODA::Profile2D> Profile2DPtr;
  typedef std::shared_ptr<YODA::Scatter1D> Scatter1DPtr;
  typedef std::shared_ptr<YODA::Scatter2D> Scatter2DPtr;
  typedef std::shared_ptr<YODA::Scatter3D> Scatter3DPtr;

}


namespace Rivet {


  /// @defgroup aotuples Minimal objects representing AO fills, to be buffered before pushToPersistent().
  ///
  /// @note Every object listed here needs a virtual fill method in YODA,
  /// otherwise the Tuple fakery won't work.
  ///
  /// @{

  /// Typedef for weights.
  using Weight = double;

  /// A single fill is a (FillType, Weight) pair.
  template <class T>
  using Fill = pair<typename T::FillType, Weight>;

  /// A set of several fill objects.
  /// @todo Why a set rather than a vector? Efficiency???
  template <class T>
  using Fills = multiset<Fill<T>>;



  /// @brief Wrappers for analysis objects to store all fills unaggregated, until collapsed by pushToPersistent().
  ///
  /// The specialisations of this inherit from the YODA analysis object types,
  /// and are used as such. The user-facing analysis objects in
  /// Analysis::analyze() are TupleWrappers on the apparent type (accessed transparently
  /// via the dereferencing of the current Wrapper<T>::active() pointer).
  ///
  /// @todo RENAME TO SOMETHING BETTER: AOProxy or FillProxy or SubEventProxy?
  template <class T>
  class TupleWrapper;


  /// TupleWrapper specialisation for Counter
  template <>
  class TupleWrapper<YODA::Counter> : public YODA::Counter {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Counter>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Counter& h) : YODA::Counter(h) {}

    /// Overloaded fill method, which stores subevent fill info until Wrapper<T>::pushToPersistent() is called.
    ///
    /// @todo Do we need to deal with users using fractions directly?
    void fill(double weight=1.0, double fraction=1.0) {
      (void)fraction; //< ???
      _fills.insert( { YODA::Counter::FillType(), weight } );
    }

    /// Empty the subevent stack (for start of new event group).
    void reset() { _fills.clear(); }

    /// Access the fill info subevent stack.
    const Fills<YODA::Counter>& fills() const { return _fills; }

  private:

    Fills<YODA::Counter> _fills;

  };


  /// TupleWrapper specialisation for Histo1D
  template <>
  class TupleWrapper<YODA::Histo1D> : public YODA::Histo1D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Histo1D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Histo1D& h) : YODA::Histo1D(h) {}

    /// Overloaded fill method, which stores subevent fill info until Wrapper<T>::pushToPersistent() is called.
    ///
    /// @todo Do we need to deal with users using fractions directly?
    void fill( double x, double weight=1.0, double fraction=1.0 ) {
      (void)fraction; //< ???
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN"); //< efficient?
      _fills.insert( { x, weight } );
    }

    /// Empty the subevent stack (for start of new event group).
    void reset() { _fills.clear(); }

    /// Access the fill info subevent stack.
    const Fills<YODA::Histo1D>& fills() const { return _fills; }

  private:

    Fills<YODA::Histo1D> _fills;

  };


  /// TupleWrapper specialisation for Profile1D
  template <>
  class TupleWrapper<YODA::Profile1D> : public YODA::Profile1D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Profile1D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Profile1D& h) : YODA::Profile1D(h) {}

    /// Overloaded fill method, which stores subevent fill info until Wrapper<T>::pushToPersistent() is called.
    ///
    /// @todo Do we need to deal with users using fractions directly?
    void fill( double x, double y, double weight=1.0, double fraction=1.0 ) {
      (void)fraction; //< ???
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN"); //< efficient?
      if ( std::isnan(y) ) throw YODA::RangeError("Y is NaN"); //< efficient?
      _fills.insert( { YODA::Profile1D::FillType{x,y}, weight } );
    }

    /// Empty the subevent stack (for start of new event group).
    void reset() { _fills.clear(); }

    /// Access the fill info subevent stack.
    const Fills<YODA::Profile1D>& fills() const { return _fills; }

  private:

    Fills<YODA::Profile1D> _fills;

  };


  /// TupleWrapper specialisation for Histo2D
  template <>
  class TupleWrapper<YODA::Histo2D> : public YODA::Histo2D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Histo2D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Histo2D& h) : YODA::Histo2D(h) {}

    /// Overloaded fill method, which stores subevent fill info until Wrapper<T>::pushToPersistent() is called.
    ///
    /// @todo Do we need to deal with users using fractions directly?
    void fill( double x, double y, double weight=1.0, double fraction=1.0 ) {
      (void)fraction; //< ???
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN"); //< efficient?
      if ( std::isnan(y) ) throw YODA::RangeError("Y is NaN"); //< efficient?
      _fills.insert( { YODA::Histo2D::FillType{x,y}, weight } );
    }

    /// Empty the subevent stack (for start of new event group).
    void reset() { _fills.clear(); }

    /// Access the fill info subevent stack.
    const Fills<YODA::Histo2D>& fills() const { return _fills; }

  private:

    Fills<YODA::Histo2D> _fills;

  };


  /// TupleWrapper specialisation for Profile2D
  template <>
  class TupleWrapper<YODA::Profile2D> : public YODA::Profile2D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Profile2D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Profile2D& h) : YODA::Profile2D(h) {}

    /// Overloaded fill method, which stores subevent fill info until Wrapper<T>::pushToPersistent() is called.
    ///
    /// @todo Do we need to deal with users using fractions directly?
    void fill( double x, double y, double z, double weight=1.0, double fraction=1.0 ) {
      (void)fraction; //< ???
      if ( std::isnan(x) ) throw YODA::RangeError("X is NaN"); //< efficient?
      if ( std::isnan(y) ) throw YODA::RangeError("Y is NaN"); //< efficient?
      if ( std::isnan(z) ) throw YODA::RangeError("Z is NaN"); //< efficient?
      _fills.insert( { YODA::Profile2D::FillType{x,y,z}, weight } );
    }

    /// Empty the subevent stack (for start of new event group).
    void reset() { _fills.clear(); }

    /// Access the fill info subevent stack.
    const Fills<YODA::Profile2D>& fills() const { return _fills; }

  private:

    Fills<YODA::Profile2D> _fills;

  };


  /// TupleWrapper specialisation for Scatter1D
  template <>
  class TupleWrapper<YODA::Scatter1D> : public YODA::Scatter1D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Scatter1D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Scatter1D& h) : YODA::Scatter1D(h) {}

  };


  /// TupleWrapper specialisation for Scatter2D
  template <>
  class TupleWrapper<YODA::Scatter2D> : public YODA::Scatter2D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Scatter2D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Scatter2D& h) : YODA::Scatter2D(h) {}

  };


  /// TupleWrapper specialisation for Scatter3D
  template <>
  class TupleWrapper<YODA::Scatter3D> : public YODA::Scatter3D {
  public:

    /// @todo Can we remove this, now that we're not relying on the AO type having a Ptr property?
    typedef shared_ptr<TupleWrapper<YODA::Scatter3D>> Ptr;

    /// @todo Can we reduce the expense of calling the full base class constructor, which mostly won't be used?
    TupleWrapper(const YODA::Scatter3D& h) : YODA::Scatter3D(h) {}

  };

  /// @}





  /// @brief Abstract interface to a set of YODA AnalysisObjects.
  ///
  /// This layer of interface is separated from the next to allow unweighted handling of
  /// Scatter objects... but do we want that? Revisit when Scatters -> Binned<Meas> and
  /// the live/dead finalize() treatment is ready to go.
  ///
  /// @todo RENAME TO SOMETHING BETTER! This is an e.g. MultiweightAOWrapper.
  ///
  /// @note This interface is not used anywhere other than in this file: eliminate/merge?
  class AnalysisObjectWrapper {
  public:

    virtual ~AnalysisObjectWrapper() {}

    /// Access the active analysis object for function calls.
    virtual YODA::AnalysisObject* operator -> () = 0;
    /// Access the active analysis object for const function calls.
    virtual YODA::AnalysisObject* operator -> () const = 0;
    /// Access the active analysis object as a reference.
    virtual const YODA::AnalysisObject& operator * () const = 0;

    /// Set active object for analyze.
    virtual void setActiveWeightIdx(size_t iWeight) = 0;

    /// Set active object for finalize.
    virtual void setActiveFinalWeightIdx(size_t iWeight) = 0;

    /// Unset the active-object pointer.
    ///
    /// @note This is for development only: we shouldn't need this in real runs.
    virtual void unsetActiveWeight() = 0;

    /// Test for equality.
    bool operator == (const AnalysisObjectWrapper& p) { return (this == &p); }
    /// Test for inequality.
    bool operator != (const AnalysisObjectWrapper& p) { return (this != &p); }

  };



  /// @todo
  /// implement scatter1dptr and scatter2dptr here
  /// these need to be multi-weighted eventually.
  /*
    class Scatter1DPtr : public AnalysisObjectPtr {
    public:
    Scatter1DPtr() : _persistent() { }

    Scatter1DPtr(size_t len_of_weightvec, const YODA::Scatter1D& p) {
    for (size_t m = 0; m < len_of_weightvec; ++m)
    _persistent.push_back(make_shared<YODA::Scatter1D>(p));
    }

    bool operator!() const { return !_persistent; }
    explicit operator bool() const { return bool(_persistent); }

    YODA::Scatter1D* operator->() { return _persistent.get(); }

    YODA::Scatter1D* operator->() const { return _persistent.get(); }

    YODA::Scatter1D & operator*() { return *_persistent; }

    const YODA::Scatter1D & operator*() const { return *_persistent; }

    protected:
    vector<YODA::Scatter1DPtr> _persistent;
    };

    class Scatter2DPtr : public AnalysisObjectPtr {
    public:
    Scatter2DPtr(size_t len_of_weightvec, const YODA::Scatter2D& p) {
    for (size_t m = 0; m < len_of_weightvec; ++m)
    _persistent.push_back(make_shared<YODA::Scatter2D>(p));
    }

    Scatter2DPtr() : _persistent() { }

    bool operator!() { return !_persistent; }
    explicit operator bool() { return bool(_persistent); }

    YODA::Scatter2D* operator->() { return _persistent.get(); }

    YODA::Scatter2D* operator->() const { return _persistent.get(); }

    YODA::Scatter2D & operator*() { return *_persistent; }

    const YODA::Scatter2D & operator*() const { return *_persistent; }

    protected:
    vector<YODA::Scatter2DPtr> _persistent;
    };

    class Scatter3DPtr : public AnalysisObjectPtr {
    public:
    Scatter3DPtr(size_t len_of_weightvec, const YODA::Scatter3D& p) {
    for (size_t m = 0; m < len_of_weightvec; ++m)
    _persistent.push_back(make_shared<YODA::Scatter3D>(p));
    }

    Scatter3DPtr() : _persistent() { }

    bool operator!() { return !_persistent; }
    explicit operator bool() { return bool(_persistent); }

    YODA::Scatter3D* operator->() { return _persistent.get(); }

    YODA::Scatter3D* operator->() const { return _persistent.get(); }

    YODA::Scatter3D & operator*() { return *_persistent; }

    const YODA::Scatter3D & operator*() const { return *_persistent; }

    protected:
    vector<YODA::Scatter3DPtr> _persistent;
    };
  */


  /// Extended abstract interface to a set of YODA AOs corresponding to multiple weight-streams, with subevent handling.
  ///
  /// @todo RENAME TO SOMETHING BETTER! This really adds the subevent proxying, so e.g. SubEventAOWrapper or MultiEventAOWrapper
  ///
  /// @note This interface is not used anywhere other than in this file: eliminate/merge?
  class MultiweightAOWrapper : public AnalysisObjectWrapper {
  public:

    /// The type being represented is a generic AO.
    using Inner = YODA::AnalysisObject;

    /// Add a new layer of subevent fill staging.
    virtual void newSubEvent() = 0;

    /// Sync the fill proxies to the persistent histogram.
    virtual void pushToPersistent(const vector<std::valarray<double> >& weight, double nlowfrac=0.0) = 0;

    /// Sync the persistent histograms to the final collection.
    virtual void pushToFinal() = 0;

    /// @todo Rename to active()?
    virtual YODA::AnalysisObjectPtr activeYODAPtr() const = 0;

    /// The histogram path, without a variation suffix.
    virtual string basePath() const = 0;

  };



  /// Type-specific multi-weight YODA analysis object wrapper.
  ///
  /// Specialisations of this (to each type of YODA object) are effectively the
  /// user-facing types in Rivet analyses, modulo a further wrapping via
  /// the Rivet shared-pointer type. They can expose either TupleWrapper<T>
  /// or T active pointers, for the analyze() and finalize() steps respectively.
  ///
  /// @todo RENAME TO SOMETHING BETTER: Wrapper<T> is far too generic. Even
  /// AnalysisObjectWrapper (with renamed base classes) would be a better
  /// user-facing choice.
  ///
  /// @todo Some things are not really well-defined here. For instance: fill()
  /// in the finalize() method and integral() in the analyze() method.
  template <class T>
  class Wrapper : public MultiweightAOWrapper {
  public:

    friend class Analysis;
    friend class AnalysisHandler;

    /// Typedef for the YODA type being represented
    using Inner = T;
    /// Typedef for a polymorphic pointer to T
    /// @note Used either to point to a real T, or a TupleWrapper<T> derived class
    using TPtr = shared_ptr<T>;


    Wrapper() = default;

    Wrapper(const vector<string>& weightnames, const T& p);

    ~Wrapper();


    /// Get the current active analysis object (may be either persistent or final, depending on stage)
    shared_ptr<T> active() const;

    /// Get the AO path of the object, without variation suffix
    string basePath() const { return _basePath; }

    /// Get the AO name of the object, without variation suffix
    string baseName() const { return _baseName; }


    /// Test for object validity.
    explicit operator bool() const { return static_cast<bool>(_active); } // Don't use active() here, assert will catch

    /// Test for object invalidity.
    bool operator ! () const { return !_active; } // Don't use active() here, assert will catch


    /// Forwarding dereference-call operator.
    T* operator -> () { return active().get(); }

    /// Forwarding dereference-call operator.
    T* operator -> () const { return active().get(); }

    /// Forwarding dereference operator.
    T& operator * () { return *active(); }

    /// Forwarding dereference operator.
    const T& operator * () const { return *active(); }


    /// Equality operator
    /// @todo These probably need to loop over all? Do we even want to provide equality? How about... no
    friend bool operator == (Wrapper a, Wrapper b){
      if (a._persistent.size() != b._persistent.size())
        return false;

      for (size_t i = 0; i < a._persistent.size(); i++) {
        if (a._persistent.at(i) != b._persistent.at(i)) {
          return false;
        }
      }
      return true;
    }

    /// Inequality operator
    friend bool operator != (Wrapper a, Wrapper b) {
      return !(a == b);
    }

    /// Less-than operator
    friend bool operator < (Wrapper a, Wrapper b) {
      if (a._persistent.size() >= b._persistent.size())
        return false;
      for (size_t i = 0; i < a._persistent.size(); i++) {
        if (*(a._persistent.at(i)) >= *(b._persistent.at(i))) {
          return false;
        }
      }
      return true;
    }


    /// Direct access to the YODA type in weight stream @a iWeight
    ///
    /// @note This is naturally a private member accessible only to , but is exposed publicly to
    /// allow analyses that explicitly study weight distributions,
    /// e.g. MC_WEIGHTS. The "private style" leading underscore in the name
    /// highlights that this should not normally be called by users.
    ///
    /// @todo Rename to minimize the clash with persistent()... or just expose persistent() instead?
    T* _getPersistent(size_t iWeight) { return _persistent.at(iWeight).get(); }



  private:

    /// @name Restricted-access methods
    ///
    /// These methods are accessible via the base classes, in which they are
    /// public, but they cannot be called directly on a Wrapper<T> object except
    /// by the friend classes Analysis (but not its children) and AnalysisHandler.
    ///
    /// @todo Review this design: it's counterintuitive and hard to maintain:
    /// how crucial is it to hide these functions from analysis authors? Can we
    /// just make them public and remove the friend stuff, then be able to use
    /// the access control properly for *really* private things?
    ///
    /// @{

    /// Set the active-object pointer to point at a variation in the persistent set
    void setActiveWeightIdx(size_t iWeight) {
      _active = _persistent.at(iWeight);
    }

    /// Set the active-object pointer to point at a variation in the final set
    void setActiveFinalWeightIdx(size_t iWeight) {
      _active = _final.at(iWeight);
    }

    /// Unset the active-object pointer
    void unsetActiveWeight() { _active.reset(); }

    /// Clear the active object pointer
    void reset() { active()->reset(); }


    /// @brief Create new object analysis-object wrappers for this sub-event
    ///
    /// Called every sub-event by AnalysisHandler::analyze() before dispatch to Analysis::analyze().
    /// The fill values will be redistributed over variations by pushToPersistent().
    void newSubEvent();

    /// Collapse the _evgroup set of tuple wrappers into fills of the persistent objects, using fractional fills if there are subevents
    void pushToPersistent(const vector<std::valarray<double> >& weight, double nlowfrac=0.0);

    /// Copy all variations from the "live" persistent set to the final collection used by Analysis::finalize()
    void pushToFinal();


    /// Get the set of persistent (i.e. after whole event groups) live objects, as used by Analysis::analyze()
    const vector<shared_ptr<T>>& persistent() const { return _persistent; }

    /// Get the set of final analysis objects, as used by Analysis::finalize() and written out
    const vector<shared_ptr<T>>& final() const { return _final; }

    /// Get the currently active (tuple wrapper on) analysis object
    virtual YODA::AnalysisObjectPtr activeYODAPtr() const { return _active; }

    /// @todo Do we need an implicit cast?
    // operator typename TPtr () { return _active; }

    /// @}


    /// @name Data members
    /// @{

    /// M of these, one for each weight
    vector<shared_ptr<T>> _persistent;

    /// The copy of M-entry _persistent that will be passed to finalize().
    vector<shared_ptr<T>> _final;

    /// A set of M subevent-wrappers, one for each weight, each containing one entry for each of the N events in evgroup.
    vector<shared_ptr<TupleWrapper<T>>> _evgroup;

    /// The currently active AO (or AO proxy).
    shared_ptr<T> _active;

    /// The base AO path of this object, without any variation suffix.
    string _basePath;

    /// The base AO name, without any variation suffix.
    string _baseName;

  };



  /// Shared-pointer type for multi-weighted Rivet AOs, dispatching through two layers of indirection.
  ///
  /// We need our own shared_ptr class, so we can dispatch -> and *
  /// all the way down to the inner YODA analysis objects
  ///
  /// @todo Provide remaining functionality that shared_ptr has (not needed right now).
  ///
  /// @todo RENAME TO SOMETHING BETTER! This naming is too generic, and the
  /// "rivet" is redundant: we need something like ao_shared_ptr or AOWrapPtr.
  template <typename T>
  class rivet_shared_ptr {
  public:
    typedef T value_type;

    rivet_shared_ptr() = default;

    rivet_shared_ptr(decltype(nullptr)) : _p(nullptr) {}

    /// Convenience constructor, pass through to the Wrapper constructor
    rivet_shared_ptr(const vector<string>& weightNames, const typename T::Inner& p)
      : _p( make_shared<T>(weightNames, p) )
    {}

    /// @todo Use SFINAE to require T<-U? Why not require rvalue == T?
    template <typename U>
    rivet_shared_ptr(const shared_ptr<U>& p)
      : _p(p)
    {}

    /// @todo Use SFINAE to require T<-U? Why not require rvalue == T?
    template <typename U>
    rivet_shared_ptr(const rivet_shared_ptr<U>& p)
      : _p(p.get())
    {}

    /// Goes right through to the active Wrapper<YODA> object's members
    T& operator -> () {
      if (_p == nullptr) throw Error("Dereferencing null AnalysisObject pointer. Is there an unbooked histogram variable?");
      return *_p;
    }

    /// Goes right through to the active Wrapper<YODA> object's members
    const T& operator -> () const                {
      if (_p == nullptr) throw Error("Dereferencing null AnalysisObject pointer. Is there an unbooked histogram variable?");
      return *_p;
    }

    /// The active YODA object
    typename T::Inner & operator * ()             { return **_p; }
    const typename T::Inner & operator * () const { return **_p; }

    /// Object validity check.
    explicit operator bool()  const { return _p && bool(*_p); }

    /// Object invalidity check.
    bool operator ! () const { return !_p || !(*_p);   }

    /// Object validity check.
    template <typename U>
    bool operator == (const rivet_shared_ptr<U>& other) const {
      return _p == other._p;
    }

    /// Object invalidity check.
    template <typename U>
    bool operator != (const rivet_shared_ptr<U>& other) const {
      return _p != other._p;
    }

    /// Less-than for ptr ordering.
    template <typename U>
    bool operator < (const rivet_shared_ptr<U>& other) const {
      return _p < other._p;
    }

    /// Greater-than for ptr ordering.
    template <typename U>
    bool operator > (const rivet_shared_ptr<U>& other) const {
      return _p > other._p;
    }

    /// Less-equals for ptr ordering.
    template <typename U>
    bool operator <= (const rivet_shared_ptr<U> & other) const {
      return _p <= other._p;
    }

    /// Greater-equals for ptr ordering.
    template <typename U>
    bool operator >= (const rivet_shared_ptr<U> & other) const {
      return _p >= other._p;
    }

    /// Get the internal shared ptr.
    shared_ptr<T> get() const { return _p; }

  private:

    /// The type being wrapped.
    shared_ptr<T> _p;

  };


  /// @defgroup useraos User-facing analysis object wrappers
  ///
  /// @note Every object listed here needs a virtual fill() method in YODA,
  /// otherwise the Tuple fakery won't work.
  ///
  /// @{

  using MultiweightAOPtr = rivet_shared_ptr<MultiweightAOWrapper>;

  using Histo1DPtr   = rivet_shared_ptr<Wrapper<YODA::Histo1D>>;
  using Histo2DPtr   = rivet_shared_ptr<Wrapper<YODA::Histo2D>>;
  using Profile1DPtr = rivet_shared_ptr<Wrapper<YODA::Profile1D>>;
  using Profile2DPtr = rivet_shared_ptr<Wrapper<YODA::Profile2D>>;
  using CounterPtr   = rivet_shared_ptr<Wrapper<YODA::Counter>>;
  using Scatter1DPtr = rivet_shared_ptr<Wrapper<YODA::Scatter1D>>;
  using Scatter2DPtr = rivet_shared_ptr<Wrapper<YODA::Scatter2D>>;
  using Scatter3DPtr = rivet_shared_ptr<Wrapper<YODA::Scatter3D>>;

  using YODA::Counter;
  using YODA::Histo1D;
  using YODA::HistoBin1D;
  using YODA::Histo2D;
  using YODA::HistoBin2D;
  using YODA::Profile1D;
  using YODA::ProfileBin1D;
  using YODA::Profile2D;
  using YODA::ProfileBin2D;
  using YODA::Scatter1D;
  using YODA::Point1D;
  using YODA::Scatter2D;
  using YODA::Point2D;
  using YODA::Scatter3D;
  using YODA::Point3D;

  ///@}





  /// @defgroup aomanip Analysis object manipulation functions
  /// @{

  /// Function to get a map of all the refdata in a paper with the
  /// given @a papername.
  map<string, YODA::AnalysisObjectPtr> getRefData(const string& papername);

  /// @todo Also provide a Scatter3D getRefData() version?

  /// Get the file system path to the reference file for this paper.
  string getDatafilePath(const string& papername);


  /// Traits class to access the type of the AnalysisObject in the reference files.
  template<typename T> struct ReferenceTraits {};
  template <> struct ReferenceTraits<Counter> { typedef Counter RefT; };
  template <> struct ReferenceTraits<Scatter1D> { typedef Scatter1D RefT; };
  template <> struct ReferenceTraits<Histo1D> { typedef Scatter2D RefT; };
  template <> struct ReferenceTraits<Profile1D> { typedef Scatter2D RefT; };
  template <> struct ReferenceTraits<Scatter2D> { typedef Scatter2D RefT; };
  template <> struct ReferenceTraits<Histo2D> { typedef Scatter3D RefT; };
  template <> struct ReferenceTraits<Profile2D> { typedef Scatter3D RefT; };
  template <> struct ReferenceTraits<Scatter3D> { typedef Scatter3D RefT; };

  /// If @a dst and @a src both are of same subclass T, copy the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false.
  template <typename T>
  inline bool aocopy(YODA::AnalysisObjectPtr src, YODA::AnalysisObjectPtr dst) {
    shared_ptr<T> tsrc = dynamic_pointer_cast<T>(src);
    if ( !tsrc ) return false;
    shared_ptr<T> tdst = dynamic_pointer_cast<T>(dst);
    if ( !tdst ) return false;
    *tdst = *tsrc;
    return true;
  }

  /// If @a dst and @a src both are of same subclass T, copy the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false. The @a scale argument will be ued to scale the weights of
  /// non-scatter types, cf. aoadd().
  template <typename T>
  inline bool aocopy(YODA::AnalysisObjectPtr src, YODA::AnalysisObjectPtr dst, double scale) {
    if (!aocopy<T>(src, dst)) return false;
    dynamic_pointer_cast<T>(dst)->scaleW(scale);
    return true;
  }

  /// If @a dst and @a src both are of same subclass T, add the
  /// contents of @a src into @a dst and return true. Otherwise return
  /// false.
  template <typename T>
  inline bool aoadd(YODA::AnalysisObjectPtr dst, YODA::AnalysisObjectPtr src, double scale) {
    shared_ptr<T> tsrc = dynamic_pointer_cast<T>(src);
    if ( !tsrc ) return false;
    shared_ptr<T> tdst = dynamic_pointer_cast<T>(dst);
    if ( !tdst ) return false;
    tsrc->scaleW(scale); //< note semi-accidental modification of the input
    try {
      *tdst += *tsrc;
    } catch (YODA::LogicError&) {
      return false;
    }
    return true;
  }

  /// If @a dst is the same subclass as @a src, copy the contents of @a
  /// src into @a dst and return true. Otherwise return false.
  bool copyao(YODA::AnalysisObjectPtr src, YODA::AnalysisObjectPtr dst, double scale=1.0);

  /// If @a dst is the same subclass as @a src, scale the contents of
  /// @a src with @a scale and add it to @a dst and return true. Otherwise
  /// return false.
  bool addaos(YODA::AnalysisObjectPtr dst, YODA::AnalysisObjectPtr src, double scale);

  /// Check if two analysis objects have the same binning or, if not
  /// binned, are in other ways compatible.
  template <typename TPtr>
  inline bool bookingCompatible(TPtr a, TPtr b) {
    return a->sameBinning(*b);
  }
  inline bool bookingCompatible(CounterPtr, CounterPtr) {
    return true;
  }
  inline bool bookingCompatible(Scatter1DPtr a, Scatter1DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(Scatter2DPtr a, Scatter2DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(Scatter3DPtr a, Scatter3DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(YODA::CounterPtr, YODA::CounterPtr) {
    return true;
  }
  inline bool bookingCompatible(YODA::Scatter1DPtr a, YODA::Scatter1DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(YODA::Scatter2DPtr a, YODA::Scatter2DPtr b) {
    return a->numPoints() == b->numPoints();
  }
  inline bool bookingCompatible(YODA::Scatter3DPtr a, YODA::Scatter3DPtr b) {
    return a->numPoints() == b->numPoints();
  }

  /// @}



  /// Class representing a YODA path with all its components.
  class AOPath {
  public:

    /// Constructor
    AOPath(string fullpath)
      : _valid(false), _path(fullpath), _raw(false), _tmp(false), _ref(false) {
      _valid = init(fullpath);
    }

    /// The full path.
    string path() const { return _path; }

    /// The analysis name.
    string analysis() const { return _analysis; }

    /// The analysis name with options.
    string analysisWithOptions() const { return _analysis + _optionstring; }

    /// The base name of the analysis object.
    string name() const { return _name; }

    /// The weight name.
    string weight() const { return _weight; }

    /// The weight component of the path
    string weightComponent() const { 
      if (_weight == "")  return _weight;
      return "[" + _weight + "]";
    }

    /// Is This a RAW (filling) object?
    bool   isRaw() const { return _raw; }

    // Is This a temporary (filling) object?
    bool   isTmp() const { return _tmp; }

    /// Is This a reference object?
    bool   isRef() const { return _ref; }

    /// The string describing the options passed to the analysis.
    string optionString() const { return _optionstring; }

    /// Are there options passed to the analysis?
    bool   hasOptions() const { return !_options.empty(); }

    /// Don't pass This optionto the analysis
    void   removeOption(string opt) { _options.erase(opt); fixOptionString(); }

    /// Pass this option to the analysis.
    void   setOption(string opt, string val) { _options[opt] = val; fixOptionString();}

    /// Was This option passed to the analyisi.
    bool   hasOption(string opt) const { return _options.find(opt) != _options.end(); }

    /// Get the value of this option.
    string getOption(string opt) const {
      auto it = _options.find(opt);
      if ( it != _options.end() ) return it->second;
      return "";
    }

    /// Reset the option string after changes;
    void fixOptionString();

    /// Creat a full path (and set) for this.
    string mkPath() const;
    string setPath() { return _path = mkPath(); }

    /// Print out information
    void debug() const;

    /// Make this class ordered.
    bool operator<(const AOPath & other) const {
      return _path < other._path;
    }

    /// Check if path is valid.
    bool valid() const { return _valid; };
    bool operator!() const { return !valid(); }

  private:

    /// Internal functions for disassembling a path name
    bool init(string fullpath);
    bool chopweight(string & fullpath);
    bool chopoptions(string & anal);

    bool _valid;
    string _path;
    string _analysis;
    string _optionstring;
    string _name;
    string _weight;
    bool _raw;
    bool _tmp;
    bool _ref;
    map<string,string> _options;

  };

}

#endif
