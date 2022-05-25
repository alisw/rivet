// -*- C++ -*-
#ifndef RIVET_PERCENTILEPROJECTION_HH
#define RIVET_PERCENTILEPROJECTION_HH

#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <map>

namespace Rivet {

/**
   @brief class for projections that reports the percentile for a
   given SingleValueProjection when initialized with a Histo1D of the
   distribution in the SingleValueProjection.

   @author Leif Lönnblad

*/
class PercentileProjection : public SingleValueProjection {
public:

  /// Constructor taking a SingleValueProjection and a calibration
  /// histogram. If increasing it means that low values corresponds to
  /// lower percentiles.
  ///
  /// @todo Use mkScatter to pass this to the Scatter2D-calibrated version?
  PercentileProjection(const SingleValueProjection & sv, const Histo1D& calhist,
                  bool increasing = false)
    : _calhist("EMPTY"), _increasing(increasing) {
    setName("PercentileProjection");
    declare(sv, "OBSERVABLE");
    //if ( !calhist ) return;
    MSG_DEBUG("Constructing PercentileProjection from " << calhist.path());
    _calhist = calhist.path();
    int N = calhist.numBins();
    double sum = calhist.sumW();

    if ( increasing ) {
      double acc = calhist.underflow().sumW();
      _table.insert(make_pair(calhist.bin(0).xEdges().first, 100.0*acc/sum));
      for ( int i = 0; i < N; ++i ) {
        acc += calhist.bin(i).sumW();
        _table.insert(make_pair(calhist.bin(i).xEdges().second, 100.0*acc/sum));
      }
    } else {
      double acc = calhist.overflow().sumW();
      _table.insert(make_pair(calhist.bin(N - 1).xEdges().second, 100.0*acc/sum));
      for ( int i = N - 1; i >= 0; --i ) {
        acc += calhist.bin(i).sumW();
        _table.insert(make_pair(calhist.bin(i).xEdges().first, 100.0*acc/sum));
      }
    }
    if (getLog().isActive(Log::DEBUG)) {
      MSG_DEBUG("Mapping from observable to percentile:");
      for (auto p : _table) {
	std::cout << std::setw(16) << p.first << " -> "
		  << std::setw(16) << p.second << "%" << std::endl;
	if (not increasing and p.second <= 0) break;
	if (increasing and p.second >= 100) break;
      }
    }
  }

  // Constructor taking a SingleValueProjection and a calibration
  // histogram. If increasing it means that low values corresponds to
  // lower percentiles.
  PercentileProjection(const SingleValueProjection & sv, const Scatter2D& calscat,
                  bool increasing = false)
    : _calhist("EMPTY"), _increasing(increasing) {
    declare(sv, "OBSERVABLE");

    //if ( !calscat ) return;
    MSG_DEBUG("Constructing PercentileProjection from " << calscat.path());
    _calhist = calscat.path();
    int N = calscat.numPoints();
    double sum = 0.0;
    for ( const auto & p : calscat.points() ) sum += p.y();

    double acc = 0.0;
    if ( increasing ) {
      _table.insert(make_pair(calscat.point(0).xMin(), 100.0*acc/sum));
      for ( int i = 0; i < N; ++i ) {
        acc += calscat.point(i).y();
        _table.insert(make_pair(calscat.point(i).xMax(), 100.0*acc/sum));
      }
    } else {
      _table.insert(make_pair(calscat.point(N - 1).xMax(), 100.0*acc/sum));
      for ( int i = N - 1; i >= 0; --i ) {
        acc += calscat.point(i).y();
        _table.insert(make_pair(calscat.point(i).xMin(), 100.0*acc/sum));
      }
    }
  }

  DEFAULT_RIVET_PROJ_CLONE(PercentileProjection);

  // The projection function takes the assigned SingeValueProjection
  // and sets the value of this projection to the corresponding
  // percentile. If no calibration has been provided, -1 will be
  // returned. If values are outside of the calibration histogram, 0
  // or 100 will be returned.
  void project(const Event& e) {
    clear();
    if ( _table.empty() ) return;
    auto&  pobs = apply<SingleValueProjection>(e, "OBSERVABLE");
    double obs  = pobs();
    double pcnt = lookup(obs);
    if ( pcnt >= 0.0 ) set(pcnt);
    MSG_DEBUG("Observable(" << pobs.name() << ")="
	      << std::setw(16) << obs 
	      << "-> Percentile=" << std::setw(16) << pcnt << "%");
  }


  // Standard comparison function.
  CmpState compare(const Projection& p) const {
    const PercentileProjection pp = dynamic_cast<const PercentileProjection&>(p);
    return mkNamedPCmp(p, "OBSERVABLE") ||
      cmp(_increasing, pp._increasing) ||
      cmp(_calhist, pp._calhist);
  }

private:

  // The (interpolated) lookup table
  double lookup(double obs) const {
    auto low = _table.upper_bound(obs);
    if ( low == _table.end() ) return _increasing? 100.0: 0.0;
    if ( low == _table.begin() ) return _increasing? 0.0: 100.0;
    auto high = low--;
    return low->second + (obs - low->first)*(high->second - low->second)/
      (high->first - low->first);
  }

  // Astring identifying the calibration histogram.
  string _calhist;

  // A lookup table to find (by interpolation) the percentile given
  // the value of the underlying SingleValueProjection.
  map<double,double> _table;

  // A flag to say whether the distribution should be integrated from
  // below or above.
  bool _increasing;

};




}

#endif
