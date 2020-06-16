// -*- C++ -*-
#ifndef RIVET_SINGLEVALUEPROJECTION_HH
#define RIVET_SINGLEVALUEPROJECTION_HH

#include "Rivet/Projection.hh"

namespace Rivet {


  /// @brief Base class for projections returning a single floating point value.
  ///
  /// @author Leif Lönnblad
  ///
  /// Project an event down to a single floating point value accessible
  /// through the operator() function.
  ///
  class SingleValueProjection: public Projection {
  public:

    /// The default constructor.
    SingleValueProjection() : _value(-1.0), _isSet(false) {
      setName("SingleValueProjection");
    }

    /// Returns true if the value has been set.
    bool isValueSet() const {
      return _isSet;
    }
    /// @deprecated Less clear alias
    bool isSet() const { return isValueSet(); }

    /// Return the single value.
    double value() const {
      return _value;
    }

    /// Return the single value.
    double operator()() const {
      return value();
    }


  protected:

    /// Set the value.
    void setValue(double v) {
      _value = v;
      _isSet = true;
    }
    /// @deprecated Less clear alias
    void set(double v) { setValue(v); }

    /// Unset the value.
    void clear() {
      _value = -1.0;
      _isSet = false;
    }


  private:

    double _value;

    /// @todo Use std::optional?
    bool _isSet;

  };


}

#endif
