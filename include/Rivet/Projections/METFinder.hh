// -*- C++ -*-
#ifndef RIVET_METFinder_HH
#define RIVET_METFinder_HH

#include "Rivet/Projection.hh"

namespace Rivet {


  /// Interface for projections that find missing transverse energy/momentum
  class METFinder : public Projection {
  public:

    /// @name Transverse momentum functions
    ///
    /// @note This may be what you want, even if the paper calls it "missing Et"!
    ///@{

    /// The vector-summed visible transverse momentum in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing pT vector.
    /// @todo Currently equivalent to vectorEt
    virtual const Vector3& vectorPt() const = 0;

    /// Convenience vector MPT function
    const Vector3 vectorMissingPt() const { return -vectorPt(); }
    // Alias
    const Vector3 vectorMPT() const { return vectorMissingPt(); }

    /// The vector-summed missing transverse momentum in the event.
    double missingPt() const { return vectorPt().mod(); }

    ///@}


    /// @name Transverse energy functions
    ///
    /// @warning Despite the common names "MET" and "SET", what's often meant is the pT functions above!
    ///@{

    /// The vector-summed visible transverse energy in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing ET vector.
    virtual const Vector3& vectorEt() const = 0;

    /// Convenience vector MET function
    const Vector3 vectorMissingEt() const { return -vectorEt(); }
    // Alias
    const Vector3 vectorMET() const { return vectorMissingEt(); }

    /// The vector-summed missing transverse energy in the event.
    double missingEt() const { return vectorEt().mod(); }
    /// Alias for missingEt
    double met() const { return missingEt(); }

    ///@}


    /// Reset the projection. Smearing functions will be unchanged.
    virtual void reset() {  }

  };


}

#endif
