// -*- C++ -*-
#ifndef RIVET_UndressBeamLeptons_HH
#define RIVET_UndressBeamLeptons_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project out the incoming beams, but subtract any colinear
  /// photons from lepton beams within a given cone.
  class UndressBeamLeptons : public Beam {
  public:

    /// Default (and only) constructor. Takes an angle as
    /// argument. The momentum of any photon within This angle wrt. a
    /// charged lepton beam will be subtracted from the beam lepton
    /// momentum.
    UndressBeamLeptons(double theta = 0.0): _thetamax(theta) {
      setName("UndressBeamLeptons");
      addProjection(FinalState(), "FS");
    }

    /// Clone on the heap
    DEFAULT_RIVET_PROJ_CLONE(UndressBeamLeptons);


    /// Project on to the Event
    virtual void project(const Event& e);


  private:

    /// Compare with other projections.
    virtual int compare(const Projection & p) const;

    /// The beam particles in the current collision
    double _thetamax;

  };


}

#endif
