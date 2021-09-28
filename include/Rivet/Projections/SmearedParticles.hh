// -*- C++ -*-
#ifndef RIVET_SmearedParticles_HH
#define RIVET_SmearedParticles_HH

#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/ParticleFinder.hh"
#include "Rivet/Tools/SmearingFunctions.hh"

namespace Rivet {


  // Recursive variadic template arg decoding
  namespace {
    template<typename T>
    vector<ParticleEffSmearFn>& toEffSmearFns(vector<ParticleEffSmearFn>& v, const T& t) {
      v.push_back(ParticleEffSmearFn(t));
      return v;
    }
    template<typename T, typename... ARGS>
    vector<ParticleEffSmearFn>& toEffSmearFns(vector<ParticleEffSmearFn>& v, const T& first, ARGS... args) {
      v.push_back(ParticleEffSmearFn(first));
      toEffSmearFns(v, args...);
      return v;
    }
  }



  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedParticles : public ParticleFinder {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with const efficiency
    SmearedParticles(const ParticleFinder& pf,
                     double eff,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {{eff}}, c)
    {    }

    /// @brief Constructor with an efficiency function
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleEffFn& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {{effFn}}, c)
    {    }

    /// @brief Constructor with const efficiency followed by a smearing function
    SmearedParticles(const ParticleFinder& pf,
                     double eff, const ParticleSmearFn& smearFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {eff, smearFn}, c)
    {    }

    /// @brief Constructor with a smearing function followed by const efficiency
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleSmearFn& smearFn, double eff,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {smearFn, eff}, c)
    {    }

    /// @brief Constructor with an efficiency function followed by a smearing function
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleEffFn& effFn, const ParticleSmearFn& smearFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {effFn, smearFn}, c)
    {    }

    /// @brief Constructor with a smearing function followed by an efficiency function
    SmearedParticles(const ParticleFinder& pf,
                     const ParticleSmearFn& smearFn, const ParticleEffFn& effFn,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, {smearFn, effFn}, c)
    {    }

    /// @brief Constructor with an ordered list of efficiency and/or smearing functions
    SmearedParticles(const ParticleFinder& pf,
                     const vector<ParticleEffSmearFn>& effSmearFns,
                     const Cut& c=Cuts::open())
      : ParticleFinder(c),
        _detFns(effSmearFns)
    {
      setName("SmearedParticles");
      declare(pf, "TruthParticles");
    }

    /// @brief Constructor with an ordered list of efficiency and/or smearing functions
    SmearedParticles(const ParticleFinder& pf,
                     const initializer_list<ParticleEffSmearFn>& effSmearFns,
                     const Cut& c=Cuts::open())
      : SmearedParticles(pf, vector<ParticleEffSmearFn>{effSmearFns}, c)
    {    }

    /// @brief Constructor with a variadic ordered list of efficiency and smearing function args
    /// @note The Cut must be provided *before* the eff/smearing functions
    /// @todo Wouldn't it be nice if the Cut could also go *after* the parameter pack?
    template<typename... ARGS>
    SmearedParticles(const ParticleFinder& pf, const Cut& c, ARGS... effSmearFns)
      : SmearedParticles(pf, toEffSmearFns(_detFns, effSmearFns...), c)
    {    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedParticles);

    //@}


    /// Compare to another SmearedParticles
    ///
    /// @note Comparing detector functions doesn't work for functors/lambdas,
    /// hence are always treated as not equivalent
    CmpState compare(const Projection& p) const {
      const SmearedParticles& other = dynamic_cast<const SmearedParticles&>(p);

      // Compare truth particles definitions
      const CmpState teq = mkPCmp(other, "TruthParticles");
      if (teq != CmpState::EQ) return teq;

      // Compare lists of detector functions
      const CmpState nfeq = cmp(_detFns.size(), other._detFns.size());
      MSG_TRACE("Numbers of detector functions = " << _detFns.size() << " VS " << other._detFns.size());
      if (nfeq != CmpState::EQ) return nfeq;
      for (size_t i = 0; i < _detFns.size(); ++i) {
        const CmpState feq = _detFns[i].cmp(other._detFns[i]);
        if (feq != CmpState::EQ) return feq;
      }

      // If we got this far, we're equal
      MSG_DEBUG("Equivalent detected! " << p.name() << ", " << this->name());
      return CmpState::EQ;
    }


    /// Perform the particle finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Particles& truthparticles = apply<ParticleFinder>(e, "TruthParticles").particlesByPt(); //truthParticles();
      _theParticles.clear(); _theParticles.reserve(truthparticles.size());
      for (const Particle& p : truthparticles) {
        Particle pdet = p;
        double peff = -1;
        bool keep = true;
        MSG_TRACE("Number of detector functions = " << _detFns.size());
        for (const ParticleEffSmearFn& fn : _detFns) {
          std::tie(pdet, peff) = fn(pdet); // smear & eff
          // Test the short-circuit random numbers if possible; note handling of < 0 and > 1 probabilities
          if (peff <= 0 || rand01() > peff) keep = false;
          MSG_DEBUG("New det particle: pid=" << pdet.pid()
                    << ", mom=" << pdet.mom()/GeV << " GeV, "
                    << "pT=" << pdet.pT()/GeV << ", eta=" << pdet.eta()
                    << " : eff=" << 100*peff << "%, discarded=" << std::boolalpha << !keep);
          if (!keep) break; // discarded; no need to try more smear-eff functions
        }
        // If discarding, go straight to the next particle
        if (!keep) continue;
        // Store, recording where the smearing was built from
        pdet.addConstituent(p); ///< @todo Is this a good idea?? What if raw particles are requested?
        _theParticles.push_back(pdet);
      }
    }

    /// Get the truth particles (sorted by pT)
    const Particles truthParticles() const {
      return getProjection<ParticleFinder>("TruthParticles").particlesByPt();
    }

    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _theParticles.clear(); }


  private:

    /// Stored efficiency & smearing functions
    vector<ParticleEffSmearFn> _detFns;

  };


}

#endif
