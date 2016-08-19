// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedParticles : public ParticleFinder {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with efficiency and smearing function args
    template <typename P2DFN>
    SmearedParticles(const ParticleFinder& pf,
                     const P2DFN& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, effFn, PARTICLE_SMEAR_IDENTITY, c)
    {    }


    /// @brief Constructor with efficiency and smearing function args
    template <typename P2DFN, typename P2PFN>
    SmearedParticles(const ParticleFinder& pf,
                     const P2DFN& effFn, const P2PFN& smearFn,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _effFn(effFn), _smearFn(smearFn)
    {
      setName("SmearedParticles");
      addProjection(pf, "TruthParticles");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedParticles);

    //@}


    /// Compare to another SmearedParticles
    int compare(const Projection& p) const {
      const SmearedParticles& other = dynamic_cast<const SmearedParticles&>(p);
      if (get_address(_effFn) == 0) return UNDEFINED;
      if (get_address(_smearFn) == 0) return UNDEFINED;
      MSG_TRACE("Eff hashes = " << get_address(_effFn) << "," << get_address(other._effFn) << "; " <<
                "smear hashes = " << get_address(_smearFn) << "," << get_address(other._smearFn));
      return mkPCmp(other, "TruthParticles") ||
        cmp(get_address(_effFn), get_address(other._effFn)) ||
        cmp(get_address(_smearFn), get_address(other._smearFn));
    }


    /// Perform the particle finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Particles& truthparticles = applyProjection<ParticleFinder>(e, "TruthParticles").particlesByPt();
      _theParticles.clear(); _theParticles.reserve(truthparticles.size());
      for (const Particle& p : truthparticles) {
        const double peff = _effFn ? _effFn(p) : 1;
        MSG_DEBUG("Efficiency of particle with pid=" << p.pid()
                  << ", mom=" << p.mom()/GeV << " GeV, "
                  << "pT=" << p.pT()/GeV << ", eta=" << p.eta()
                  << " : " << 100*peff << "%");
        if (peff <= 0) continue; //< no need to roll expensive dice (and we deal with -ve probabilities, just in case)
        if (peff < 1 && rand01() > peff) continue; //< roll dice (and deal with >1 probabilities, just in case)
        _theParticles.push_back(_smearFn ? _smearFn(p) : p); //< smearing
      }
    }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _theParticles.clear(); }


  private:

    /// Stored efficiency function
    std::function<double(const Particle&)> _effFn;

    /// Stored smearing function
    std::function<Particle(const Particle&)> _smearFn;

  };


}

#endif
