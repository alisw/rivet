// -*- C++ -*-
#ifndef RIVET_ConstLossyFinalState_HH
#define RIVET_ConstLossyFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Random.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LossyFinalState.hh"

namespace Rivet {


  /// Functor used to implement constant random lossiness.
  class ConstRandomFilter {
  public:

    ConstRandomFilter(double lossFraction)
      : _lossFraction(lossFraction)
    {
      assert(_lossFraction >= 0);
    }

    // If operator() returns true, particle is deleted ("lost")
    bool operator()(const Particle&) {
      return rand01() < _lossFraction;
    }

    CmpState compare(const ConstRandomFilter& other) const {
      return cmp(_lossFraction, other._lossFraction);
    }

  private:

    double _lossFraction;

  };



  /// @brief Randomly lose a constant fraction of particles.
  class ConstLossyFinalState : public LossyFinalState<ConstRandomFilter> {
  public:

    /// @name Constructors
    //@{

    /// Constructor from a FinalState.
    ConstLossyFinalState(const FinalState& fsp, double lossfraction)
      : LossyFinalState<ConstRandomFilter>(fsp, ConstRandomFilter(lossfraction))
    {
      setName("ConstLossyFinalState");
    }

    /// Stand-alone constructor. Initialises the base FinalState projection.
    ConstLossyFinalState(double lossfraction, const Cut& c=Cuts::open())
      : LossyFinalState<ConstRandomFilter>(ConstRandomFilter(lossfraction), c)
    {
      setName("ConstLossyFinalState");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(ConstLossyFinalState);

    //@}

  };


}

#endif
