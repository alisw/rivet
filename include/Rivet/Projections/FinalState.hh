// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Project out all final-state particles in an event.
  /// Probably the most important projection in Rivet!
  class FinalState : public ParticleFinder {
  public:

    /// @name Standard constructors etc.
    //@{

    /// Construction using Cuts object
    FinalState(const Cut& c=Cuts::open());

    /// Construction using another FinalState and a Cuts object
    FinalState(const FinalState& fsp, const Cut& c);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(FinalState);

    //@}


    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// Compare projections.
    virtual CmpState compare(const Projection& p) const;

    /// Decide if a particle is to be accepted or not.
    /// @todo Rename to _accept or acceptFinal?
    virtual bool accept(const Particle& p) const;


  private:

    // Hide lossy copy constructors for all classes derived from FinalState
    template<typename T> FinalState(const T& rhs);
    template<typename T> FinalState const& operator=(T const& rhs);

  };


}

#endif
