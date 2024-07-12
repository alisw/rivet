// -*- C++ -*-
#ifndef RIVET_GammaGammaFinalState_HH
#define RIVET_GammaGammaFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/GammaGammaKinematics.hh"

namespace Rivet {


  /// @brief Final state particles boosted to the hadronic center of mass system.
  ///
  /// NB. The GammaGamma scattered leptons are not included in the final state particles.
  class GammaGammaFinalState: public FinalState {
  public:

    /// @name Constructors
    //@{
    /// Constructor with optional FinalState
    /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    GammaGammaFinalState(const FinalState& fs=FinalState(), const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
    {
      setName("GammaGammaFinalState");
      declare(fs, "FS");
      declare(kinematicsp, "Kinematics");
    }

    /// Constructor with explicit cuts to define final-state particles
    /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    GammaGammaFinalState(const Cut& c, const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
      : GammaGammaFinalState(FinalState(c), kinematicsp)
    {    }

    // /// @brief Constructor with default FinalState
    // /// @note The GammaGammaKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    // GammaGammaFinalState(const GammaGammaKinematics& kinematicsp=GammaGammaKinematics())
    //   : GammaGammaFinalState(FinalState(), kinematicsp)
    // {    }

    /// Backward compatible constructor with default FinalState
    /// @deprecated Prefer a version that doesn't need a GammaGammaKinematics argument
    GammaGammaFinalState(const GammaGammaKinematics& kinematicsp)
      : GammaGammaFinalState(FinalState(), kinematicsp)
    {    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(GammaGammaFinalState);

    //@}


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const {
      return mkNamedPCmp(p, "Kinematics") || mkNamedPCmp(p, "FS");
    }


  };


}

#endif
