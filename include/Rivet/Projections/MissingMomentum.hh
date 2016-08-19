// -*- C++ -*-
#ifndef RIVET_MissingMomentum_HH
#define RIVET_MissingMomentum_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Calculate missing \f$ E \f$, \f$ E_\perp \f$ etc.
  ///
  /// Project out the total visible energy vector, allowing missing
  /// \f$ E \f$, \f$ E_\perp \f$ etc. to be calculated. Final state
  /// visibility restrictions are automatic.
  class MissingMomentum : public Projection {
  public:

    /// Default constructor with optional cut.
    MissingMomentum(const Cut& c=Cuts::open()) {
      setName("MissingMomentum");
      FinalState fs(c);
      addProjection(fs, "FS");
      addProjection(VisibleFinalState(fs), "VisibleFS");
    }


    /// Constructor.
    MissingMomentum(const FinalState& fs) {
      setName("MissingMomentum");
      addProjection(fs, "FS");
      addProjection(VisibleFinalState(fs), "VisibleFS");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(MissingMomentum);


    /// The vector-summed visible four-momentum in the event.
    ///
    /// @note Reverse this vector with .reverse() to get the missing momentum vector.
    ///
    /// @note The optional @a mass argument is used to set a mass on the 4-vector. By
    ///   default it is zero (since missing momentum is really a 3-momentum quantity:
    ///   adding the E components of visible momenta just gives a huge mass)
    const FourMomentum visibleMomentum(double mass=0*GeV) const;
    /// Alias for visibleMomentum
    const FourMomentum visibleMom(double mass=0*GeV) const { return visibleMomentum(mass); }

    /// The missing four-momentum in the event, required to balance the final state.
    ///
    /// @note The optional @a mass argument is used to set a mass on the 4-vector. By
    ///   default it is zero (since missing momentum is really a 3-momentum quantity:
    ///   adding the E components of visible momenta just gives a huge mass)
    const FourMomentum missingMomentum(double mass=0*GeV) const { return visibleMomentum(mass).reverse(); }
    /// Alias for missingMomentum
    const FourMomentum missingMom(double mass=0*GeV) const { return missingMomentum(mass); }

    /// The vector-summed visible transverse energy in the event, as a 3-vector with z=0
    /// @note Reverse this vector with operator- to get the missing ET vector.
    const Vector3& vectorEt() const { return _vet; }

    /// The vector-summed missing transverse energy in the event.
    double missingEt() const { return vectorEt().mod(); }
    /// Alias for missingEt
    double met() const { return missingEt(); }

    /// The scalar-summed visible transverse energy in the event.
    double scalarEt() const { return _set; }
    /// Alias for scalarEt
    double set() const { return scalarEt(); }


  protected:

    /// Apply the projection to the event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  public:

    /// Clear the projection results.
    void clear();


  private:

    /// The total visible momentum
    FourMomentum _momentum;

    /// Scalar transverse energy
    double _set;

    /// Vector transverse energy
    Vector3 _vet;

  };


}

#endif
