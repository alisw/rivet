#ifndef RIVET_PARTICLEUTILS_HH
#define RIVET_PARTICLEUTILS_HH

#include "Rivet/Particle.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

// Macros to map Rivet::Particle functions to PID:: functions of the same name
/// @todo Can leave return type out of the macro and put that on each line where it's used?
#define PARTICLE_TO_PID_BOOLFN(fname) inline bool fname (const Particle& p) { return PID:: fname (p.pid()); }
#define PARTICLE_TO_PID_INTFN(fname) inline int fname (const Particle& p) { return PID:: fname (p.pid()); }
#define PARTICLE_TO_PID_DBLFN(fname) inline double fname (const Particle& p) { return PID:: fname (p.pid()); }

namespace Rivet {


  /// @name Particle classifier functions
  //@{

  /// Unbound function access to PID code
  inline int pid(const Particle& p) { return p.pid(); }

  /// Unbound function access to abs PID code
  inline int abspid(const Particle& p) { return p.abspid(); }


  /// Is this particle species charged?
  PARTICLE_TO_PID_BOOLFN(isCharged)

  /// Is this particle species neutral?
  PARTICLE_TO_PID_BOOLFN(isNeutral)


  /// Is this a neutrino?
  PARTICLE_TO_PID_BOOLFN(isNeutrino)

  /// Determine if the PID is that of a charged lepton
  PARTICLE_TO_PID_BOOLFN(isChLepton)

  /// Determine if the PID is that of a photon
  PARTICLE_TO_PID_BOOLFN(isPhoton)

  /// Determine if the PID is that of an electron or positron
  PARTICLE_TO_PID_BOOLFN(isElectron)

  /// Determine if the PID is that of an muon or antimuon
  PARTICLE_TO_PID_BOOLFN(isMuon)

  /// Determine if the PID is that of an tau or antitau
  PARTICLE_TO_PID_BOOLFN(isTau)

  /// Determine if the PID is that of a hadron
  PARTICLE_TO_PID_BOOLFN(isHadron)

  /// Determine if the PID is that of a meson
  PARTICLE_TO_PID_BOOLFN(isMeson)

  /// Determine if the PID is that of a baryon
  PARTICLE_TO_PID_BOOLFN(isBaryon)

  /// Determine if the PID is that of a quark
  PARTICLE_TO_PID_BOOLFN(isQuark)

  /// Determine if the PID is that of a parton (quark or gluon)
  PARTICLE_TO_PID_BOOLFN(isParton)



  /// Determine if the PID is that of a W+
  PARTICLE_TO_PID_BOOLFN(isWplus)

  /// Determine if the PID is that of a W-
  PARTICLE_TO_PID_BOOLFN(isWminus)

  /// Determine if the PID is that of a W+-
  PARTICLE_TO_PID_BOOLFN(isW)

  /// Determine if the PID is that of a Z0
  PARTICLE_TO_PID_BOOLFN(isZ)

  /// Determine if the PID is that of an SM/lightest SUSY Higgs
  PARTICLE_TO_PID_BOOLFN(isHiggs)

  /// Determine if the PID is that of an s/sbar
  PARTICLE_TO_PID_BOOLFN(isStrange)

  /// Determine if the PID is that of a c/cbar
  PARTICLE_TO_PID_BOOLFN(isCharm)

  /// Determine if the PID is that of a b/bbar
  PARTICLE_TO_PID_BOOLFN(isBottom)

  /// Determine if the PID is that of a t/tbar
  PARTICLE_TO_PID_BOOLFN(isTop)


  /// Determine if the particle is a heavy flavour hadron or parton
  PARTICLE_TO_PID_BOOLFN(isHeavyFlavour)

  /// Determine if the PID is that of a heavy parton (c,b,t)
  PARTICLE_TO_PID_BOOLFN(isHeavyParton)

  /// Determine if the PID is that of a light parton (u,d,s)
  PARTICLE_TO_PID_BOOLFN(isLightParton)


  /// Determine if the PID is that of a heavy flavour (b or c) meson
  PARTICLE_TO_PID_BOOLFN(isHeavyMeson)

  /// Determine if the PID is that of a heavy flavour (b or c) baryon
  PARTICLE_TO_PID_BOOLFN(isHeavyBaryon)

  /// Determine if the PID is that of a heavy flavour (b or c) hadron
  PARTICLE_TO_PID_BOOLFN(isHeavyHadron)


  /// Determine if the PID is that of a light flavour (not b or c) meson
  PARTICLE_TO_PID_BOOLFN(isLightMeson)

  /// Determine if the PID is that of a light flavour (not b or c) baryon
  PARTICLE_TO_PID_BOOLFN(isLightBaryon)

  /// Determine if the PID is that of a light flavour (not b or c) hadron
  PARTICLE_TO_PID_BOOLFN(isLightHadron)


  /// Determine if the PID is that of a b-meson.
  PARTICLE_TO_PID_BOOLFN(isBottomMeson)

  /// Determine if the PID is that of a b-baryon.
  PARTICLE_TO_PID_BOOLFN(isBottomBaryon)

  /// Determine if the PID is that of a b-hadron.
  PARTICLE_TO_PID_BOOLFN(isBottomHadron)


  /// @brief Determine if the PID is that of a c-meson.
  ///
  /// Specifically, the _heaviest_ quark is a c: a B_c is a b-meson and NOT a c-meson.
  /// Charmonia (closed charm) are counted as c-mesons here.
  PARTICLE_TO_PID_BOOLFN(isCharmMeson)

  /// @brief Determine if the PID is that of a c-baryon.
  ///
  /// Specifically, the _heaviest_ quark is a c: a baryon containing a b & c
  /// is a b-baryon and NOT a c-baryon. To test for the simpler case, just use
  /// a combination of hasCharm() and isBaryon().
  PARTICLE_TO_PID_BOOLFN(isCharmBaryon)

  /// Determine if the PID is that of a c-hadron.
  PARTICLE_TO_PID_BOOLFN(isCharmHadron)


  // /// Determine if the PID is that of a strange meson
  // PARTICLE_TO_PID_BOOLFN(isStrangeMeson)

  // /// Determine if the PID is that of a strange baryon
  // PARTICLE_TO_PID_BOOLFN(isStrangeBaryon)

  // /// Determine if the PID is that of a strange hadron
  // PARTICLE_TO_PID_BOOLFN(isStrangeHadron)



  /// Is this a pomeron, odderon, or generic reggeon?
  PARTICLE_TO_PID_BOOLFN(isReggeon)

  /// Determine if the PID is that of a diquark (used in hadronization models)
  PARTICLE_TO_PID_BOOLFN(isDiquark)

  /// Determine if the PID is that of a pentaquark (hypothetical hadron)
  PARTICLE_TO_PID_BOOLFN(isPentaquark)

  /// Is this a fundamental SUSY particle?
  PARTICLE_TO_PID_BOOLFN(isSUSY)

  /// Is this an R-hadron?
  PARTICLE_TO_PID_BOOLFN(isRhadron)

  /// Is this a technicolor particle?
  PARTICLE_TO_PID_BOOLFN(isTechnicolor)

  /// Is this an excited (composite) quark or lepton?
  PARTICLE_TO_PID_BOOLFN(isExcited)

  /// Is this a Kaluza-Klein excitation?
  PARTICLE_TO_PID_BOOLFN(isKK)

  /// Is this a graviton?
  PARTICLE_TO_PID_BOOLFN(isGraviton)

  /// Is this a BSM particle (including graviton)?
  PARTICLE_TO_PID_BOOLFN(isBSM)



  /// Determine if the PID is in the generator-specific range
  PARTICLE_TO_PID_BOOLFN(isGenSpecific)

  /// Determine if the PID is that of an EW scale resonance
  PARTICLE_TO_PID_BOOLFN(isResonance)

  /// Check the PID for usability in transport codes like Geant4
  PARTICLE_TO_PID_BOOLFN(isTransportable)



  /// Does this particle contain an up quark?
  PARTICLE_TO_PID_BOOLFN(hasUp)

  /// Does this particle contain a down quark?
  PARTICLE_TO_PID_BOOLFN(hasDown)

  /// Does this particle contain a strange quark?
  PARTICLE_TO_PID_BOOLFN(hasStrange)

  /// Does this particle contain a charm quark?
  PARTICLE_TO_PID_BOOLFN(hasCharm)

  /// Does this particle contain a bottom quark?
  PARTICLE_TO_PID_BOOLFN(hasBottom)

  /// Does this particle contain a top quark?
  PARTICLE_TO_PID_BOOLFN(hasTop)



  /// jSpin returns 2J+1, where J is the total spin
  PARTICLE_TO_PID_INTFN(jSpin)

  /// sSpin returns 2S+1, where S is the spin
  PARTICLE_TO_PID_INTFN(sSpin)

  /// lSpin returns 2L+1, where L is the orbital angular momentum
  PARTICLE_TO_PID_INTFN(lSpin)


  /// Return the charge
  PARTICLE_TO_PID_DBLFN(charge)

  /// Return 3 times the charge (3 x quark charge is an int)
  PARTICLE_TO_PID_INTFN(charge3)

  /// Return the absolute charge
  PARTICLE_TO_PID_DBLFN(abscharge)

  /// Return 3 times the abs charge (3 x quark charge is an int)
  PARTICLE_TO_PID_INTFN(abscharge3)

  /// Alias for charge3
  /// @deprecated Use charge3
  PARTICLE_TO_PID_INTFN(threeCharge)


  /// Get the atomic number (number of protons) in a nucleus/ion
  PARTICLE_TO_PID_INTFN(nuclZ)

  /// Get the atomic weight (number of nucleons) in a nucleus/ion
  PARTICLE_TO_PID_INTFN(nuclA)

  /// If this is a nucleus (ion), get nLambda
  PARTICLE_TO_PID_INTFN(nuclNlambda)

  //@}


  /// @name Particle charge/sign comparison functions
  //@{

  /// @brief Return true if Particles @a a and @a b have the opposite charge sign
  /// @note Two neutrals returns false
  inline bool oppSign(const Particle& a, const Particle& b) {
    return sign(a.charge3()) == -sign(b.charge3()) && sign(a.charge3()) != ZERO;
  }

  /// Return true if Particles @a a and @a b have the same charge sign
  /// @note Two neutrals returns true
  inline bool sameSign(const Particle& a, const Particle& b) {
    return sign(a.charge3()) == sign(b.charge3());
  }

  /// Return true if Particles @a a and @a b have the exactly opposite charge
  /// @note Two neutrals returns false
  inline bool oppCharge(const Particle& a, const Particle& b) {
    return a.charge3() == -b.charge3() && a.charge3() != 0;
  }

  /// Return true if Particles @a a and @a b have the same charge (including neutral)
  /// @note Two neutrals returns true
  inline bool sameCharge(const Particle& a, const Particle& b) {
    return a.charge3() == b.charge3();
  }

  /// Return true if Particles @a a and @a b have a different (not necessarily opposite) charge
  inline bool diffCharge(const Particle& a, const Particle& b) {
    return a.charge3() != b.charge3();
  }

  //@}


  /// @name Particle classifying functors
  ///
  /// To be passed to any() or all() e.g. any(p.children(), HasPID(PID::MUON))
  //@{

  /// Base type for Particle -> bool functors
  struct BoolParticleFunctor {
    virtual bool operator()(const Particle& p) const = 0;
  };

  /// PID matching functor
  struct hasPID : public BoolParticleFunctor {
    hasPID(PdgId pid) : targetpid(pid) { }
    bool operator()(const Particle& p) const { return p.pid() == targetpid; }
    PdgId targetpid;
  };
  using HasPID = hasPID;

  /// |PID| matching functor
  struct hasAbsPID : public BoolParticleFunctor {
    hasAbsPID(PdgId pid) : targetpid(abs(pid)) { }
    bool operator()(const Particle& p) const { return p.abspid() == abs(targetpid); }
    PdgId targetpid;
  };
  using HasAbsPID = hasAbsPID;

  //@}



  //////////////////////////////////////



  /// @name Non-PID particle properties, via unbound functions
  //@{

  /// Is this particle potentially visible in a detector?
  inline bool isVisible(const Particle& p) { return p.isVisible(); }

  /// Decide if a given particle is prompt, via Particle::isPrompt()
  inline bool isPrompt(const Particle& p, bool from_prompt_tau=false, bool from_prompt_mu=false) {
    return p.isPrompt(from_prompt_tau, from_prompt_mu);
  }

  /// Decide if a given particle is stable, via Particle::isStable()
  inline bool isStable(const Particle& p) { return p.isStable(); }

  /// Check whether a given PID is found in the particle's ancestor list
  inline bool hasAncestor(const Particle& p, PdgId pid)  { return p.hasAncestor(pid); }

  /// Determine whether the particle is from a b-hadron decay
  inline bool fromBottom(const Particle& p) { return p.fromBottom(); }

  /// @brief Determine whether the particle is from a c-hadron decay
  inline bool fromCharm(const Particle& p) { return p.fromCharm(); }

  /// @brief Determine whether the particle is from a hadron decay
  inline bool fromHadron(const Particle& p) { return p.fromHadron(); }

  /// @brief Determine whether the particle is from a tau decay
  inline bool fromTau(const Particle& p, bool prompt_taus_only=false) {
    return p.fromTau(prompt_taus_only);
  }

  /// @brief Determine whether the particle is from a prompt tau decay
  inline bool fromPromptTau(const Particle& p) { return p.fromPromptTau(); }

  /// @brief Determine whether the particle is from a hadron or tau decay
  inline bool fromDecay(const Particle& p) { return p.fromDecay(); }

  //@}



  /// @name Unbound functions for filtering particles
  //@{

  /// Get a subset of the supplied particles that passes the supplied Cut
  Particles filterBy(const Particles& particles, const Cut& c);

  /// Filter a particle collection in-place to the subset that passes the supplied Cut
  Particles& ifilterBy(Particles& particles, const Cut& c);

  /// Filter a particle collection to the subset that passes the supplied Cut, into a new container
  /// @note New container will be replaced, not appended to
  inline Particles& filterBy(Particles& particles, const Cut& c, Particles& out) {
    //const Particles& const_particles = particles;
    out = filterBy(particles, c);
    return out;
  }

  //@}



  /// @name Particle pair functions
  //@{

  /// Get the PDG ID codes of a ParticlePair
  /// @todo Make ParticlePair a custom class instead?
  inline PdgIdPair pids(const ParticlePair& pp) {
    return make_pair(pp.first.pid(), pp.second.pid());
  }

  //@}


}

#endif
