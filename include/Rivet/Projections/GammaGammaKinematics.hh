// -*- C++ -*-
#ifndef RIVET_GammaGammaKinematics_HH
#define RIVET_GammaGammaKinematics_HH

#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/GammaGammaLeptons.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Get the gamma gamma kinematic variables and relevant boosts for an event.
  class GammaGammaKinematics : public Projection {
  public:

    /// The default constructor.
    GammaGammaKinematics(const GammaGammaLeptons & lepton = GammaGammaLeptons(),
                  const std::map<std::string,std::string> & opts =
                  std::map<std::string,std::string>())
      : _theQ2(make_pair(-1.0,-1.0)), _theW2(-1.0) //,_theX(-1.0), _theY(-1.0), _theS(-1.0)
    {
      setName("GammaGammaKinematics");
      //addPdgIdPair(ANY, hadid);
      declare(Beam(), "Beam");
      declare(lepton, "Lepton");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(GammaGammaKinematics);


  protected:

    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual CmpState compare(const Projection& p) const;


  public:

    /// The \f$Q^2\f$.
    pair<double,double> Q2() const { return _theQ2; }

    /// The \f$W^2\f$.
    double W2() const { return _theW2; }

    /// The incoming lepton beam particle
    const ParticlePair& beamLeptons() const {
      return _inLepton;
    }

    /// The scattered GammaGamma lepton
    const ParticlePair & scatteredLeptons() const {
      return _outLepton;
    }



  private:

    /// The \f$Q^2\f$.
    pair<double,double> _theQ2;

    /// The \f$W^2\f$.
    double _theW2;

    /// Incoming and outgoing GammaGamma particles
    ParticlePair _inLepton, _outLepton;

  };


}

#endif
