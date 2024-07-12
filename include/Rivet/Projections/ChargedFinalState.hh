// -*- C++ -*-
#ifndef RIVET_ChargedFinalState_HH
#define RIVET_ChargedFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project only charged final state particles.
  class ChargedFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Construction from another FinalState
    ChargedFinalState(const FinalState& fsp);

    /// Construction using Cuts object
    ChargedFinalState(const Cut& c=Cuts::open());

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(ChargedFinalState);

    //@}


    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;

  };


}


#endif
