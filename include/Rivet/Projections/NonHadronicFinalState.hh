// -*- C++ -*-
#ifndef RIVET_NonHadronicFinalState_HH
#define RIVET_NonHadronicFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project only hadronic final state particles.
  class NonHadronicFinalState : public FinalState {
  public:

    /// Constructor: the supplied FinalState projection is assumed to live through the run.
    NonHadronicFinalState(FinalState& fsp)
    {
      setName("NonHadronicFinalState");
      declare(fsp, "FS");
    }

    NonHadronicFinalState(const Cut& c=Cuts::open())
    {
      setName("NonHadronicFinalState");
      declare(FinalState(c), "FS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(NonHadronicFinalState);


    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

  };


}


#endif
