// -*- C++ -*-
#ifndef RIVET_DISDiffHadron_HH
#define RIVET_DISDiffHadron_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Get the incoming and outgoing hadron in a diffractive ep
  /// event.
  class DISDiffHadron : public Projection {
  public:

    /// @name Constructors.
    //@{

    /// Default constructor.
    DISDiffHadron() {
      setName("DISDiffHadron");
      addProjection(Beam(), "Beam");
      addProjection(FinalState(), "FS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(DISDiffHadron);

    //@}


  protected:

    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const;


  public:

    /// The incoming lepton
    const Particle& in() const { return _incoming; }

    /// The outgoing lepton
    const Particle& out() const { return _outgoing; }

  private:

    /// The incoming lepton
    Particle _incoming;

    /// The outgoing lepton
    Particle _outgoing;

  };

}


#endif
