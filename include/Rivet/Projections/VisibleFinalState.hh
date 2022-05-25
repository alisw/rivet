// -*- C++ -*-
#ifndef RIVET_VisibleFinalState_HH
#define RIVET_VisibleFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Final state modifier excluding particles which are not experimentally visible
  class VisibleFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor with min and max pseudorapidity \f$ \eta \f$ and min \f$ p_T \f$ (in GeV).
    VisibleFinalState(const Cut& c=Cuts::open())
    {
      setName("VisibleFinalState");
      declare(FinalState(c), "FS");
    }

    /// Constructor with specific FinalState.
    VisibleFinalState(const FinalState& fsp)
    {
      setName("VisibleFinalState");
      declare(fsp, "FS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(VisibleFinalState);

    //@}


    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

  };


}

#endif
