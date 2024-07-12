// -*- C++ -*-
#ifndef RIVET_InitialQuarks_HH
#define RIVET_InitialQuarks_HH

#ifndef I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#warning "This is a dangerous projection for a few specific old analyses. Not for general use!"
#endif

#include "Rivet/Projection.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Project out quarks from the hard process in \f$ e^+ e^- \to Z^0 \f$ events
  ///
  /// @deprecated We're not sure exactly when we'lll get rid of this, but it's going to happen...
  ///
  /// @warning This is a very dangerous and specific projection!
  class InitialQuarks : public Projection {
  public:

    /// @name Standard constructors and destructors.
    //@{
    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    InitialQuarks() {
      setName("InitialQuarks");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(InitialQuarks);

    //@}

    /// Access the projected final-state particles.
    virtual const Particles& particles() const { return _theParticles; }

    /// Is this final state empty?
    virtual bool empty() const { return _theParticles.empty(); }


  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// Compare projections.
    virtual CmpState compare(const Projection& p) const;


  protected:

    /// The final-state particles.
    Particles _theParticles;

  };

}


#endif
