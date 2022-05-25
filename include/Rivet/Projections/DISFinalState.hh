// -*- C++ -*-
#ifndef RIVET_DISFinalState_HH
#define RIVET_DISFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Final state particles boosted to the hadronic center of mass system.
  ///
  /// NB. The DIS scattered lepton is not included in the final state particles.
  class DISFinalState: public FinalState {
  public:

    /// Type of DIS boost to apply
    enum class BoostFrame { HCM, BREIT, LAB };


    /// @name Constructors
    //@{

    /// @brief Constructor with explicit FinalState
    ///
    /// @deprecated The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(const FinalState& fs, BoostFrame boosttype, const DISKinematics& kinematicsp=DISKinematics())
      : _boosttype(boosttype)
    {
      setName("DISFinalState");
      declare(fs, "FS");
      declare(kinematicsp, "Kinematics");
    }

    /// @brief Constructor with optional FinalState
    ///
    /// @deprecated The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(BoostFrame boosttype, const FinalState& fs=FinalState(), const DISKinematics& kinematicsp=DISKinematics())
      : DISFinalState(fs, boosttype, kinematicsp)
    {    }

    /// @brief Constructor with explicit cuts to define final-state particles
    ///
    /// @note The cuts will be applied *before* the boost, e.g. to express detector acceptance.
    ///
    /// @todo Add a second optional Cut argument for post-boost cuts.
    ///
    /// @deprecated The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(const Cut& c, BoostFrame boosttype, const DISKinematics& kinematicsp=DISKinematics())
      : DISFinalState(FinalState(c), boosttype, kinematicsp)
    {    }

    /// @brief Constructor with explicit cuts to define final-state particles
    ///
    /// @note The cuts will be applied *before* the boost, e.g. to express detector acceptance.
    ///
    /// @todo Add a second optional Cut argument for post-boost cuts.
    ///
    /// @deprecated The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(BoostFrame boosttype, const Cut& c, const DISKinematics& kinematicsp=DISKinematics())
      : DISFinalState(FinalState(c), boosttype, kinematicsp)
    {    }

    // /// @brief Constructor with default FinalState
    // ///
    // /// @note The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    // DISFinalState(BoostFrame boosttype, const DISKinematics& kinematicsp=DISKinematics())
    //   : DISFinalState(FinalState(), boosttype, kinematicsp)
    // {    }

    /// @brief Backward-compatible constructor with default FinalState
    ///
    /// @deprecated Prefer a version that doesn't need a DISKinematics argument
    DISFinalState(const DISKinematics& kinematicsp, BoostFrame boosttype)
      : DISFinalState(FinalState(), boosttype, kinematicsp)
    {    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(DISFinalState);

    //@}


    /// Get the associated DISKinematics (to avoid needing a separate projection)
    const DISKinematics& kinematics() {
      return getProjection<DISKinematics>("Kinematics");
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const {
      const DISFinalState& other = dynamic_cast<const DISFinalState&>(p);
      return mkNamedPCmp(p, "Kinematics") || mkNamedPCmp(p, "FS") || cmp(_boosttype, other._boosttype);
    }


  private:

    BoostFrame _boosttype;

  };


}

#endif
