// -*- C++ -*-
#ifndef RIVET_MissingMomentum_HH
#define RIVET_MissingMomentum_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/METFinder.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Calculate missing \f$ E \f$, \f$ E_\perp \f$ etc. as complements to the total visible momentum
  ///
  /// Project out the total visible energy vector, allowing missing \f$ E \f$,
  /// \f$ E_\perp \f$ etc. to be calculated. Final-state particle visibility
  /// restrictions are automatic, and the resulting visible/missing momentum
  /// vectors are over the whole event rather than over hard objects (jets +
  /// leptons) or specific to prompt invisibles.
  class MissingMomentum : public METFinder {
  public:

    /// Canonical constructor taking a FinalState as argument
    MissingMomentum(const FinalState& fs) {
      setName("MissingMomentum");
      declare(fs, "FS");
      declare(VisibleFinalState(fs), "VisibleFS");
    }

    /// Default constructor with optional cut
    MissingMomentum(const Cut& c=Cuts::open())
      : MissingMomentum(FinalState(c))
    {    }


    /// Clone on the heap
    DEFAULT_RIVET_PROJ_CLONE(MissingMomentum);


    /// @name Visible/missing four-momentum functions
    /// @{

    /// @brief The vector-summed visible four-momentum in the event
    ///
    /// @note Reverse this vector with .reverse() to get the missing momentum vector.
    ///
    /// @note The optional @a mass argument is used to set a mass on the 4-vector. By
    ///   default it is zero (since missing momentum is really a 3-momentum quantity:
    ///   adding the E components of visible momenta just gives a huge mass)
    const FourMomentum visibleMomentum(double mass=0*GeV) const;
    /// Alias for visibleMomentum
    const FourMomentum visibleMom(double mass=0*GeV) const { return visibleMomentum(mass); }

    /// @brief The missing four-momentum in the event, required to balance the final state.
    ///
    /// @note The optional @a mass argument is used to set a mass on the 4-vector. By
    ///   default it is zero (since missing momentum is really a 3-momentum quantity:
    ///   adding the E components of visible momenta just gives a huge mass)
    const FourMomentum missingMomentum(double mass=0*GeV) const { return visibleMomentum(mass).reverse(); }
    /// Alias for missingMomentum
    const FourMomentum missingMom(double mass=0*GeV) const { return missingMomentum(mass); }

    /// @}


    /// @name Transverse momentum functions
    ///
    /// @note This may be what you want, even if the paper calls it "missing Et"!
    /// @{

    /// @brief The vector-summed visible transverse momentum in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing pT vector.
    const Vector3& vectorPt() const { return _vpt; }

    /// The scalar-summed visible transverse momentum in the event.
    double scalarPt() const { return _spt; }
    // /// Alias for scalarPt
    // double spt() const { return scalarPt(); }

    /// @}


    /// @name Transverse energy functions
    ///
    /// @warning Despite the common names "MET" and "SET", what's often meant is the pT functions above!
    /// @{

    /// @brief The vector-summed visible transverse energy in the event, as a 3-vector with z=0
    ///
    /// @note Reverse this vector with operator- to get the missing ET vector.
    const Vector3& vectorEt() const { return _vet; }

    /// The scalar-summed visible transverse energy in the event.
    double scalarEt() const { return _set; }
    /// Alias for scalarEt
    double set() const { return scalarEt(); }

    /// @}


    /// Clear the projection results.
    void clear();


  protected:

    /// Apply the projection to the event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  private:

    /// The total visible momentum
    FourMomentum _momentum;

    /// Scalar transverse energy
    double _set, _spt;

    /// Vector transverse energy
    Vector3 _vet, _vpt;

  };


}

#endif
