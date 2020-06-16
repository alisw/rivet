// -*- C++ -*-
#ifndef RIVET_CENTRALITYPROJECTION_HH
#define RIVET_CENTRALITYPROJECTION_HH

#include "Rivet/Projections/PercentileProjection.hh"
#include "Rivet/Tools/RivetYODA.hh"
#include <map>

namespace Rivet {

/// @brief Used together with the percentile-based analysis objects Percentile and PercentileXaxis
///
/// The interior actually defines several different centrality
/// estimates: the centrality observable used in the experiment with a
/// reference calibration ("REF"); the same but using a user-defined
/// calibration done with the corresponding minimum bias analysis
/// ("GEN"); a centrality based on the impact parameter reported in
/// HepMC::HeavyIon::impact_parameter, using a calibration histogram
/// generated with the same minimum bias analysis ("IMP"). For HepMC3
/// it may optionally also include a direct report from the generator
/// about the centrality, if available in HepMC::HeavyIon::centrality
/// ("RAW"), and a user-defined generated centrality estimate
/// communicated via the HepMC::HeavyIon::user_cent_estimate ("USR").
///
/// @author Leif Lönnblad
class CentralityProjection: public SingleValueProjection {
public:

  /// Default constructor
  CentralityProjection() {}


  DEFAULT_RIVET_PROJ_CLONE(CentralityProjection);


  /// @brief Add a new centrality estimate.
  ///
  /// The SingleValueProjection, @a p, should return a value between 0
  /// and 100, and the @a pname should be one of "REF", "GEN", "IMP",
  /// "USR", or "RAW", as described above.
  void add(const SingleValueProjection & p, string pname) {
    _projNames.push_back(pname);
    declare(p, pname);
  }

  /// Perform all internal projections.
  void project(const Event& e) {
    _values.clear();
    for ( string pname : _projNames )
      _values.push_back(apply<SingleValueProjection>(e, pname)());
    if ( !_values.empty() ) set(_values[0]);
  }

  /// Cheek if no internal projections have been added.
  bool empty() const {
    return _projNames.empty();
  }

  /// Return the percentile of the @a i'th projection.
  ///
  /// Note that operator() will return the zero'th projection.
  double operator[](int i) const {
    return _values[i];
  }

  // Standard comparison function.
  CmpState compare(const Projection& p) const {
    const CentralityProjection* other = dynamic_cast<const CentralityProjection*>(&p);
    if (other->_projNames.size() == 0) return CmpState::NEQ;
    for (string pname : _projNames) {
      bool hasPname = true;
      for (string p2name : other->_projNames){
        if (pname != p2name) hasPname = false;
      }
      if (!hasPname) return CmpState::NEQ;
    }
    return CmpState::EQ;
  }

  /// The list of names of the internal projections.
  vector<string> projections() const {
    return _projNames;
  }

private:

  /// The list of names of the internal projections.
  vector<string> _projNames;

  /// The list of percentiles resulting from the last projection.
  vector<double> _values;

};

}

#endif
