// -*- C++ -*-
#ifndef RIVET_SINGLEVALUEPROJECTION_HH
#define RIVET_SINGLEVALUEPROJECTION_HH

#include "Rivet/Projection.hh"

namespace Rivet {

/** @brief Base class for projections returning a single floating
    point value.

    @author Leif Lönnblad

    Project an event down to a single floating point value accessible
    through the operator() function.

*/
class SingleValueProjection: public Projection {

public:

  /// The default constructor.
  SingleValueProjection() : _value(-1.0), _isSet(false) {
    setName("SingleValueProjection");
  }

  /// Returns true if the value has been set.
  bool isSet() const {
    return _isSet;
  }

  /// Return the single value.
  double operator()() const {
    return _value;
  }

protected:

  /// Set the value.
  void set(double v) {
    _value = v;
    _isSet = true;
  }


  /// Unset the value.
  void clear() {
    _value = -1.0;
    _isSet = false;
  }

private:

  double _value;

  bool _isSet;

};

}

#endif

