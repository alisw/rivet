// -*- C++ -*-
#ifndef RIVET_LossyFinalState_HH
#define RIVET_LossyFinalState_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Templated FS projection which can lose some of the supplied particles.
  template <typename FILTER>
  class LossyFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor from FinalState.
    LossyFinalState(const FinalState& fsp, FILTER filter)
      : _filter(filter)
    {
      setName("LossyFinalState");
      declare(fsp, "FS");
    }

    /// Stand-alone constructor. Initialises the base FinalState projection.
    LossyFinalState(FILTER filter, const Cut& c=Cuts::open())
      : _filter(filter)
    {
      setName("LossyFinalState");
      declare(FinalState(c), "FS");
    }

    /// Virtual destructor, to allow subclassing
    virtual ~LossyFinalState() { }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(LossyFinalState);

    //@}


    /// Apply the projection on the supplied event.
    void project(const Event& e) {
      const FinalState& fs = applyProjection<FinalState>(e, "FS");
      getLog() << Log::DEBUG << "Pre-loss number of FS particles = " << fs.particles().size() << '\n';
      _theParticles.clear();
      std::remove_copy_if(fs.particles().begin(), fs.particles().end(),
                          std::back_inserter(_theParticles), _filter);
      getLog() << Log::DEBUG << "Filtered number of FS particles = " << _theParticles.size() << '\n';
    }


    /// Compare projections.
    CmpState compare(const Projection& p) const {
      const LossyFinalState<FILTER>& other = pcast< LossyFinalState<FILTER> >(p);
      const CmpState fscmp = mkNamedPCmp(other, "FS");
      if (fscmp != CmpState::EQ) return fscmp;
      return _filter.compare(other._filter);
    }


  protected:

    /// Filtering object: must support operator(const Particle&) and compare(const Filter&)
    FILTER _filter;

  };


}

#endif
