#ifndef PERCENTILE_HH
#define PERCENTILE_HH

#include "Rivet/Event.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/ProjectionApplier.hh"

namespace Rivet {


  /// Forward declaration.
  class Analysis;


  /// @brief PercentileBase is the base class of all Percentile classes.
  ///
  /// This base class contains all non-templated variables and
  /// infrastructure needed.
  class PercentileBase {
  public:

    /// @brief Constructor
    ///
    /// Constructor requiring a pointer, @a ana, to the Analysis to which this
    /// object belongs and the name of the CentralityProjection, @a
    /// projname, to be used.
    PercentileBase(Analysis * ana, string projName)
      : _ana(ana), _projName(projName) {}

    /// @brief Default constructor.
    PercentileBase() {}

    /// @brief Initialize the PercentileBase for a new event.
    ///
    /// This will perform the assigned CentralityProjection and select
    /// out the (indices) of the internal AnalysisObjects that are to be
    /// active in this event.
    void selectBins(const Event &);

    /// @brief Helper function to check if @a x is within @a range.
    static bool inRange(double x, pair<float,float> range) {
      return x >= range.first && ( x < range.second || ( x == 100.0 && x == range.second ) );
    }

    /// @brief Copy information from @a other PercentileBase
    void copyFrom(const PercentileBase & other) {
      _ana = other._ana;
      _projName = other._projName;
      _cent = other._cent;
    }

    /// @brief check if @a other PercentileBase is compatible with this.
    bool compatible(const PercentileBase & other) const {
      return ( _ana == other._ana &&
               _projName == other._projName &&
               _cent == other._cent );
    }

    /// @brief return the list of centrality bins.
    ///
    /// The size of this vector is the same as number of internal
    /// analysis objects in the sub class PercentileTBase.
    const vector< pair<float, float> > & centralities() const {
      return _cent;
    }


  protected:

    /// The Analysis object to which This object is assigned.
    Analysis* _ana;

    /// The name of the CentralityProjection.
    string _projName;

    /// The list of indices of the analysis objects that are to be
    /// filled in the current event.
    vector<int> _activeBins;

    /// The list of centrality intervals, one for each included analysis
    /// object.
    vector<pair<float, float> > _cent;

  };



  /// @brief PercentileTBase is the base class of all Percentile classes.
  ///
  /// This base class contains all template-dependent variables and
  /// infrastructure needed for Percentile and PercentileXaxis.
  template<class T>
  class PercentileTBase : public PercentileBase {
  public:

    /// Convenient typedef.
    typedef rivet_shared_ptr<Wrapper<T>> TPtr;

    /// @brief Main constructor
    ///
    /// requiring a pointer, @a ana, to the Analysis to which this
    /// object belongs and the name of the CentralityProjection, @a
    /// projname, to be used.
    PercentileTBase(Analysis * ana, string projName)
      : PercentileBase(ana, projName), _histos() {}

    /// @brief Default constructor
    PercentileTBase() {}

    /// @brief Empty destructor
    ~PercentileTBase() {}

    /// @brief Add a new percentile bin
    ///
    /// Add an analysis objects which are clones of @a temp that should
    /// be active for events in the given centrality bin @a
    /// cent. Several analysis objects may be added depending on the
    /// number of alternative centrality definitions in the
    /// CentralityProjection @a proj. This function is common for
    /// Percentile and PecentileXaxis, but for the latter the @a cent
    /// argument should be left to its default.
    void add(TPtr ao, CounterPtr cnt,
             pair<float,float> cent = {0.0, 100.0} ) {
      _cent.push_back(cent);
      _histos.push_back( { ao, cnt } );
    }

    /// @brief Copy the information from an @a other Percentile object.
    ///
    /// This function differs from a simple assignement as the @a other
    /// analysis objects are not copied, but supplied separately through
    /// @a tv.
    bool add(const PercentileBase & other, const vector<TPtr> & tv) {
      copyFrom(other);
      if ( tv.size() != _cent.size() ) return false;
      for ( auto t : tv )
        _histos.push_back( { t, CounterPtr() } );
      return true;
    }

    /// @brief Initialize for a new event. Select which AnalysisObjects
    /// should be filled for this event. Keeps track of the number of
    /// events seen for each centrality bin and AnalysisAbject.
    bool init(const Event & event) {
      selectBins(event);
      for (const auto bin : _activeBins)
        _histos[bin].second->fill();
      return !_activeBins.empty();
    }

    /// @brief Normalize each AnalysisObject
    ///
    /// Normalize by dividing by the sum of the events seen for each centrality
    /// bin.
    void normalizePerEvent() {
      for (const auto &hist : _histos)
        if ( hist.second->numEntries() > 0 && hist.first->numEntries() > 0)
          hist.first->scaleW(1./hist.second->val());
    }

    /// Simple scaling of each AnalysisObject
    void scale(float scale) {
      for (const auto hist : _histos)
        hist.first->scaleW(scale);
    }

    /// Execute a function for each AnalysisObject
    void exec(function<void(T&)> f) { for ( auto hist : _histos) f(hist); }

    /// @brief Access the underlyng AnalysisObjects
    ///
    /// The returned vector contains a pair, where the first member is
    /// the AnalysisObject and the second is a counter keeping track of
    /// the sum of event weights for which the AnalysisObject has been
    /// active.
    const vector<pair<TPtr, CounterPtr > > &
    analysisObjects() const{
      return _histos;
    }


  protected:

    /// The returned vector contains a pair, where the first member is
    /// the AnalysisObject and the second is a counter keeping track of
    /// the sum of event weights for which the AnalysisObject has been
    /// active.
    vector<pair<TPtr, CounterPtr > > _histos;

  };



  /// @brief The Percentile class for centrality binning.
  ///
  /// The Percentile class automatically handles the selection of which
  /// AnalysisObject(s) should be filled depending on the centrality of
  /// an event. It cointains a list of AnalysisObjects, one for each
  /// centrality bin requested (note that these bins may be overlapping)
  /// and each centrality definition is available in the assigned
  /// CentralityProjection.
  template<class T>
  class Percentile : public PercentileTBase<T> {
  public:

    /// @brief Main constructor
    ///
    /// Requires a pointer, @a ana, to the Analysis to which this
    /// object belongs and the name of the CentralityProjection, @a
    /// projname, to be used.
    Percentile(Analysis * ana, string projName)
      : PercentileTBase<T>(ana, projName) {}

    /// @brief Default constructor.
    Percentile() {}

    /// @brief Empty destructor.
    ~Percentile() {}

    /// Needed to access members of the templated base class.
    using PercentileTBase<T>::_histos;

    /// Needed to access members of the templated base class.
    using PercentileTBase<T>::_activeBins;

    /// Fill each AnalysisObject selected in the last call to
    /// PercentileTBase<T>init
    template<typename... Args>
    void fill(Args... args) {
      for (const auto bin : _activeBins) {
        _histos[bin].first->fill(args...);
      }
    }

    /// Subtract the contents fro another Pecentile.
    Percentile<T> &operator-=(const Percentile<T> &rhs) {
      const int nCent = _histos.size();
      for (int iCent = 0; iCent < nCent; ++iCent) {
        *_histos[iCent].first -= *rhs._histos[iCent].first;
      }
    }

    /// Add the contents fro another Pecentile.
    Percentile<T> &operator+=(const Percentile<T> &rhs) {
      const int nCent = _histos.size();
      for (int iCent = 0; iCent < nCent; ++iCent) {
        *_histos[iCent].first += *rhs._histos[iCent].first;
        /// @todo should this also add the Counter?
      }
    }

    /// Make this object look like a pointer.
    Percentile<T> *operator->() { return this; }

    /// Pointer to member operator.
    Percentile<T> &operator->*(function<void(T&)> f) { exec(f);  return *this; }

  };



  /// @brief The PercentileXaxis class for centrality binning.
  ///
  /// The PercentileXaxis class automatically handles the x-axis of an
  /// AnalysisObject when the x-axis is to be the centrality of an
  /// event. This could also be done by eg. filling directly a Histo1D
  /// with the result of a CentralityProjection. However, since the
  /// CentralityProjection may handle several centrality definitions at
  /// the same time it is reasonable to instead use
  /// PercentileXaxis<Histo1D> which will fill one histogram for each
  /// centrality definition.
  ///
  /// Operationally this class works like the Percentile class, but only
  /// one centrality bin (0-100) is included. When fill()ed the first
  /// argument is always given by the assigned CentralityProjection.
  template<class T>
  class PercentileXaxis : public PercentileTBase<T> {
  public:

    /// @brief Main constructor
    ///
    /// Requires a pointer, @a ana, to the Analysis to which this
    /// object belongs and the name of the CentralityProjection, @a
    /// projname, to be used.
    PercentileXaxis(Analysis * ana, string projName)
      : PercentileTBase<T>(ana, projName) {}

    /// @brief Default constructor
    PercentileXaxis() {}

    /// @brief Empty destructor
    ~PercentileXaxis() {}

    /// Needed to access members of the templated base class
    using PercentileTBase<T>::_histos;

    /// Needed to access members of the templated base class
    using PercentileTBase<T>::_activeBins;

    /// Fill each AnalysisObject selected in the last call to
    /// PercentileTBase<T>init
    template<typename... Args>
    void fill(Args... args) {
      for (const auto bin : _activeBins) {
        _histos[bin].first->fill(bin, args...);
      }
    }

    /// Subtract the contents from another PercentileXaxis.
    PercentileXaxis<T> &operator-=(const PercentileXaxis<T> &rhs) {
      const int nCent = _histos.size();
      for (int iCent = 0; iCent < nCent; ++iCent) {
        *_histos[iCent].first -= *rhs._histos[iCent].first;
      }
    }

    /// Add the contents from another PercentileXaxis.
    PercentileXaxis<T> &operator+=(const PercentileXaxis<T> &rhs) {
      const int nCent = this->_histos.size();
      for (int iCent = 0; iCent < nCent; ++iCent) {
        *_histos[iCent].first += *rhs._histos[iCent].first;
      }
    }

    /// Make this object look like a pointer.
    PercentileXaxis<T> *operator->() { return this; }

    /// Pointer to member operator.
    PercentileXaxis<T> &operator->*(function<void(T&)> f) { exec(f);  return *this; }

  };



  /// @name Combining Percentiles
  ///
  /// Follows the naming of functions for the underlying AnalysisObjects: global operators
  // @{

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  divide(const Percentile<T> numer, const Percentile<T> denom) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( numer.compatible(denom) );
    for ( int i = 0, N = numer.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(divide(*numer.analysisObjects()[i].first,
                                                   *denom.analysisObjects()[i].first)));
    ret.add(numer, scatters);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  divide(const Percentile<T> numer,
         const Percentile<typename ReferenceTraits<T>::RefT> denom) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( numer.compatible(denom) );
    for ( int i = 0, N = numer.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(divide(*numer.analysisObjects()[i].first,
                                                   *denom.analysisObjects()[i].first)));
    ret.add(numer, scatters);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  divide(const Percentile<typename ReferenceTraits<T>::RefT> numer,
         const Percentile<T> denom) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<typename ReferenceTraits<T>::RefT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( numer.compatible(denom) );
    for ( int i = 0, N = numer.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(divide(*numer.analysisObjects()[i].first,
                                                   *denom.analysisObjects()[i].first)));
    ret.add(numer, scatters);
    return ret;
  }

  template <typename T>
  Percentile<T> add(const Percentile<T> pctla, const Percentile<T> pctlb) {
    Percentile<T> ret;
    vector<typename T::Ptr> aos;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      aos.push_back(make_shared<T>(add(*pctla.analysisObjects()[i].first,
                                       *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, aos);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  add(const Percentile<T> pctla,
      const Percentile<typename ReferenceTraits<T>::RefT> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(add(*pctla.analysisObjects()[i].first,
                                                *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  add(const Percentile<typename ReferenceTraits<T>::RefT> pctla,
      const Percentile<T> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(add(*pctla.analysisObjects()[i].first,
                                                *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  Percentile<T> subtract(const Percentile<T> pctla, const Percentile<T> pctlb) {
    Percentile<T> ret;
    vector<typename T::Ptr> aos;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      aos.push_back(make_shared<T>(subtract(*pctla.analysisObjects()[i].first,
                                            *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, aos);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  subtract(const Percentile<T> pctla,
           const Percentile<typename ReferenceTraits<T>::RefT> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(subtract(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  subtract(const Percentile<typename ReferenceTraits<T>::RefT> pctla,
           const Percentile<T> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(subtract(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  multiply(const Percentile<T> pctla,
           const Percentile<typename ReferenceTraits<T>::RefT> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(multiply(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  multiply(const Percentile<typename ReferenceTraits<T>::RefT> pctla,
           const Percentile<T> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    Percentile<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(multiply(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }



  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  divide(const PercentileXaxis<T> numer, const PercentileXaxis<T> denom) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( numer.compatible(denom) );
    for ( int i = 0, N = numer.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(divide(*numer.analysisObjects()[i].first,
                                                   *denom.analysisObjects()[i].first)));
    ret.add(numer, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  divide(const PercentileXaxis<T> numer,
         const PercentileXaxis<typename ReferenceTraits<T>::RefT> denom) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( numer.compatible(denom) );
    for ( int i = 0, N = numer.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(divide(*numer.analysisObjects()[i].first,
                                                   *denom.analysisObjects()[i].first)));
    ret.add(numer, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  divide(const PercentileXaxis<typename ReferenceTraits<T>::RefT> numer,
         const PercentileXaxis<T> denom) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<typename ReferenceTraits<T>::RefT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( numer.compatible(denom) );
    for ( int i = 0, N = numer.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(divide(*numer.analysisObjects()[i].first,
                                                   *denom.analysisObjects()[i].first)));
    ret.add(numer, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<T> add(const PercentileXaxis<T> pctla, const PercentileXaxis<T> pctlb) {
    PercentileXaxis<T> ret;
    vector<typename T::Ptr> aos;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      aos.push_back(make_shared<T>(add(*pctla.analysisObjects()[i].first,
                                       *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, aos);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  add(const PercentileXaxis<T> pctla,
      const PercentileXaxis<typename ReferenceTraits<T>::RefT> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(add(*pctla.analysisObjects()[i].first,
                                                *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  add(const PercentileXaxis<typename ReferenceTraits<T>::RefT> pctla,
      const PercentileXaxis<T> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(add(*pctla.analysisObjects()[i].first,
                                                *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<T> subtract(const PercentileXaxis<T> pctla, const PercentileXaxis<T> pctlb) {
    PercentileXaxis<T> ret;
    vector<typename T::Ptr> aos;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      aos.push_back(make_shared<T>(subtract(*pctla.analysisObjects()[i].first,
                                            *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, aos);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  subtract(const PercentileXaxis<T> pctla,
           const PercentileXaxis<typename ReferenceTraits<T>::RefT> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(subtract(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  subtract(const PercentileXaxis<typename ReferenceTraits<T>::RefT> pctla,
           const PercentileXaxis<T> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(subtract(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  multiply(const PercentileXaxis<T> pctla,
           const PercentileXaxis<typename ReferenceTraits<T>::RefT> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(multiply(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  multiply(const PercentileXaxis<typename ReferenceTraits<T>::RefT> pctla,
           const PercentileXaxis<T> pctlb) {
    typedef typename ReferenceTraits<T>::RefT ScatT;
    PercentileXaxis<ScatT> ret;
    vector<typename ScatT::Ptr> scatters;
    assert( pctla.compatible(pctlb) );
    for ( int i = 0, N = pctla.analysisObjects().size(); i < N; ++i )
      scatters.push_back(make_shared<ScatT>(multiply(*pctla.analysisObjects()[i].first,
                                                     *pctlb.analysisObjects()[i].first)));
    ret.add(pctla, scatters);
    return ret;
  }

  template <typename T>
  Percentile<T>
  operator+(const Percentile<T> pctla, const Percentile<T> pctlb) {
    return add(pctla, pctlb);
  }

  template <typename T>
  Percentile<T>
  operator-(const Percentile<T> pctla, const Percentile<T> pctlb) {
    return subtract(pctla, pctlb);
  }

  template <typename T>
  Percentile<typename ReferenceTraits<T>::RefT>
  operator/(const Percentile<T> numer, const Percentile<T> denom) {
    return divide(numer, denom);
  }

  template <typename T>
  PercentileXaxis<T>
  operator+(const PercentileXaxis<T> pctla, const PercentileXaxis<T> pctlb) {
    return add(pctla, pctlb);
  }

  template <typename T>
  PercentileXaxis<T>
  operator-(const PercentileXaxis<T> pctla, const PercentileXaxis<T> pctlb) {
    return subtract(pctla, pctlb);
  }

  template <typename T>
  PercentileXaxis<typename ReferenceTraits<T>::RefT>
  operator/(const PercentileXaxis<T> numer, const PercentileXaxis<T> denom) {
    return divide(numer, denom);
  }


}

#endif
