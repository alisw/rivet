// -*- C++ -*-
#ifndef RIVET_ParticleSmearingFunctions_HH
#define RIVET_ParticleSmearingFunctions_HH

#include "Rivet/Particle.hh"
#include "Rivet/Tools/MomentumSmearingFunctions.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @ingroup smearing
  /// @{

  /// @defgroup smearing_particle Generic particle filtering, efficiency and smearing utils
  /// @{

  /// Typedef for Particle smearing functions/functors
  typedef function<Particle(const Particle&)> ParticleSmearFn;

  /// Typedef for Particle efficiency functions/functors
  typedef function<double(const Particle&)> ParticleEffFn;


  /// Take a Particle and return 0
  inline double PARTICLE_EFF_ZERO(const Particle& ) { return 0; }
  /// Alias for PARTICLE_EFF_ZERO
  inline double PARTICLE_EFF_0(const Particle& ) { return 0; }
  /// Alias for PARTICLE_EFF_ZERO
  inline double PARTICLE_FN0(const Particle& ) { return 0; }

  /// Take a Particle and return 1
  inline double PARTICLE_EFF_ONE(const Particle& ) { return 1; }
  /// Alias for PARTICLE_EFF_ONE
  inline double PARTICLE_EFF_1(const Particle& ) { return 1; }
  /// Alias for PARTICLE_EFF_ONE
  inline double PARTICLE_EFF_PERFECT(const Particle& ) { return 1; }
  /// Alias for PARTICLE_EFF_ONE
  inline double PARTICLE_FN1(const Particle& ) { return 1; }

  /// Take a Particle and return a constant number
  struct PARTICLE_EFF_CONST {
    PARTICLE_EFF_CONST(double x) : _x(x) {}
    double operator () (const Particle& )  const { return _x; }
    double _x;
  };


  /// Take a Particle and return it unmodified
  inline Particle PARTICLE_SMEAR_IDENTITY(const Particle& p) { return p; }
  /// Alias for PARTICLE_SMEAR_IDENTITY
  inline Particle PARTICLE_SMEAR_PERFECT(const Particle& p) { return p; }


  /// @brief Functor for simultaneous efficiency-filtering and smearing of Particles
  ///
  /// A central element of the SmearedParticles system
  struct ParticleEffSmearFn {
    ParticleEffSmearFn(const ParticleSmearFn& s, const ParticleEffFn& e)
      : sfn(s), efn(e) {    }

    ParticleEffSmearFn(const ParticleEffFn& e, const ParticleSmearFn& s)
      : sfn(s), efn(e) {    }

    ParticleEffSmearFn(const ParticleSmearFn& s)
      : sfn(s), efn(PARTICLE_EFF_ONE) {    }

    ParticleEffSmearFn(const ParticleEffFn& e)
      : sfn(PARTICLE_SMEAR_IDENTITY), efn(e) {    }

    ParticleEffSmearFn(double eff)
      : ParticleEffSmearFn(PARTICLE_EFF_CONST(eff)) {    }

    /// Smear and calculate an efficiency for the given particle
    pair<Particle,double> operator() (const Particle& p) const {
      return make_pair(sfn(p), efn(p));
    }

    /// Compare to another, for use in the projection system
    CmpState cmp(const ParticleEffSmearFn& other) const {
      // cout << "Eff hashes = " << get_address(efn) << "," << get_address(other.efn) << "; "
      //      << "smear hashes = " << get_address(sfn) << "," << get_address(other.sfn) << '\n';
      if (get_address(sfn) == 0 || get_address(other.sfn) == 0) return CmpState::NEQ;
      if (get_address(efn) == 0 || get_address(other.efn) == 0) return CmpState::NEQ;
      return Rivet::cmp(get_address(sfn), get_address(other.sfn)) || Rivet::cmp(get_address(efn), get_address(other.efn));
    }

    /// Automatic conversion to a smearing function
    operator ParticleSmearFn () { return sfn; }
    /// Automatic conversion to an efficiency function
    operator ParticleEffFn () { return efn; }

    // Stored functions/functors
    const ParticleSmearFn sfn;
    const ParticleEffFn efn;
  };


  /// Return true if Particle @a p is chosen to survive a random efficiency selection
  inline bool efffilt(const Particle& p, const ParticleEffFn& feff) {
    return rand01() < feff(p);
  }

  /// @brief A functor to return true if Particle @a p survives a random efficiency selection
  ///
  /// @deprecated Prefer... ?
  struct ParticleEffFilter {
    template <typename FN>
    ParticleEffFilter(const FN& feff) : _feff(feff) {}
    ParticleEffFilter(double eff) : ParticleEffFilter( [&](const Particle& p){return eff;} ) {}
    bool operator () (const Particle& p)  const { return efffilt(p, _feff); }
  private:
    const ParticleEffFn _feff;
  };
  using particleEffFilter = ParticleEffFilter;

  /// @}

  /// @}

}

#endif
