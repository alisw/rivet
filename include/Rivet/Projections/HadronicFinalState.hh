// -*- C++ -*-
#ifndef RIVET_HadronicFinalState_HH
#define RIVET_HadronicFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project only hadronic final state particles.
  class HadronicFinalState : public FinalState {
  public:

    /// Constructor: the supplied FinalState projection is assumed to live through the run.
    HadronicFinalState(const FinalState& fsp)
    {
      setName("HadronicFinalState");
      declare(fsp, "FS");
    }

    HadronicFinalState(const Cut& c=Cuts::open())
    {
      setName("HadronicFinalState");
      declare(FinalState(c), "FS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(HadronicFinalState);

  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

  };


}


#endif
