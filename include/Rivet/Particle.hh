// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Particle.fhh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Math/LorentzTrans.hh"
// NOTE: Rivet/Tools/ParticleUtils.hh included at the end
#include "fastjet/PseudoJet.hh"

namespace Rivet {


  /// Particle representation, either from a HepMC::GenEvent or reconstructed.
  class Particle : public ParticleBase {
  public:

    /// @name Constructors
    //@{

    /// Default constructor.
    /// @note A particle without info is useless. This only exists to keep STL containers happy.
    Particle()
      : ParticleBase(),
        _original(0), _id(0)
    { }

    /// Constructor without GenParticle.
    Particle(PdgId pid, const FourMomentum& mom, const FourVector& pos=FourVector())
      : ParticleBase(),
        _original(0), _id(pid), _momentum(mom), _origin(pos)
    { }

    /// Constructor from a HepMC GenParticle pointer.
    Particle(const GenParticle* gp)
      : ParticleBase(),
        _original(gp), _id(gp->pdg_id()),
        _momentum(gp->momentum())
    {
      const GenVertex* vprod = gp->production_vertex();
      if (vprod != NULL)
        setOrigin(vprod->position().t(), vprod->position().x(), vprod->position().y(), vprod->position().z());
    }

    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp)
      : ParticleBase(),
        _original(&gp), _id(gp.pdg_id()),
        _momentum(gp.momentum())
    {
      const GenVertex* vprod = gp.production_vertex();
      if (vprod != NULL)
        setOrigin(vprod->position().t(), vprod->position().x(), vprod->position().y(), vprod->position().z());
    }

    //@}


    /// @name Other representations and implicit casts
    //@{

    /// Converter to FastJet3 PseudoJet
    virtual fastjet::PseudoJet pseudojet() const {
      return fastjet::PseudoJet(mom().px(), mom().py(), mom().pz(), mom().E());
    }

    /// Cast operator to FastJet3 PseudoJet
    operator PseudoJet () const { return pseudojet(); }


    /// Get a const pointer to the original GenParticle.
    const GenParticle* genParticle() const {
      return _original;
    }

    /// Cast operator for conversion to GenParticle*
    operator const GenParticle* () const { return genParticle(); }

    //@}


    /// @name Kinematic properties
    //@{

    /// The momentum.
    const FourMomentum& momentum() const {
      return _momentum;
    }

    /// Set the momentum.
    Particle& setMomentum(const FourMomentum& momentum) {
      _momentum = momentum;
      return *this;
    }

    /// Set the momentum via components.
    Particle& setMomentum(double E, double px, double py, double pz) {
      _momentum = FourMomentum(E, px, py, pz);
      return *this;
    }

    /// Apply an active Lorentz transform to this particle
    Particle& transformBy(const LorentzTransform& lt);

    //@


    /// @name Positional properties
    //@{

    /// The origin position.
    const FourVector& origin() const {
      return _origin;
    }
    /// Set the origin position.
    Particle& setOrigin(const FourVector& position) {
      _origin = position;
      return *this;
    }
    /// Set the origin position via components.
    Particle& setOrigin(double t, double x, double y, double z) {
      _origin = FourMomentum(t, x, y, z);
      return *this;
    }

    //@}


    /// @name Particle ID code accessors
    //@{

    /// This Particle's PDG ID code.
    PdgId pid() const { return _id; }
    /// Absolute value of the PDG ID code.
    PdgId abspid() const { return std::abs(_id); }
    /// This Particle's PDG ID code (alias).
    /// @deprecated Prefer the pid/abspid form
    PdgId pdgId() const { return _id; }

    //@}


    /// @name Particle species properties
    //@{

    /// The charge of this Particle.
    double charge() const { return PID::charge(pid()); }

    /// The absolute charge of this Particle.
    double abscharge() const { return PID::abscharge(pid()); }

    /// Three times the charge of this Particle (i.e. integer multiple of smallest quark charge).
    int charge3() const { return PID::charge3(pid()); }

    /// Alias for charge3
    /// @deprecated Use charge3
    int threeCharge() const { return PID::threeCharge(pid()); }

    /// Three times the absolute charge of this Particle (i.e. integer multiple of smallest quark charge).
    int abscharge3() const { return PID::abscharge3(pid()); }

    /// Is this a hadron?
    bool isHadron() const { return PID::isHadron(pid()); }

    /// Is this a meson?
    bool isMeson() const { return PID::isMeson(pid()); }

    /// Is this a baryon?
    bool isBaryon() const { return PID::isBaryon(pid()); }

    /// Is this a lepton?
    bool isLepton() const { return PID::isLepton(pid()); }

    /// Is this a neutrino?
    bool isNeutrino() const { return PID::isNeutrino(pid()); }

    /// Does this (hadron) contain a b quark?
    bool hasBottom() const { return PID::hasBottom(pid()); }

    /// Does this (hadron) contain a c quark?
    bool hasCharm() const { return PID::hasCharm(pid()); }

    // /// Does this (hadron) contain an s quark?
    // bool hasStrange() const { return PID::hasStrange(pid()); }

    /// Is this particle potentially visible in a detector?
    bool isVisible() const;

    //@}


    /// @name Ancestry properties
    //@{

    /// Check whether a given PID is found in the particle's ancestor list
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasAncestor(PdgId pdg_id) const;

    /// @brief Determine whether the particle is from a b-hadron decay
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromBottom() const;

    /// @brief Determine whether the particle is from a c-hadron decay
    ///
    /// @note If a hadron contains b and c quarks it is considered a bottom
    /// hadron and NOT a charm hadron.
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromCharm() const;

    // /// @brief Determine whether the particle is from a s-hadron decay
    // ///
    // /// @note If a hadron contains b or c quarks as well as strange it is
    // /// considered a b or c hadron, but NOT a strange hadron.
    // ///
    // /// @note This question is valid in MC, but may not be perfectly answerable
    // /// experimentally -- use this function with care when replicating
    // /// experimental analyses!
    // bool fromStrange() const;

    /// @brief Determine whether the particle is from a hadron decay
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromHadron() const;

    /// @brief Determine whether the particle is from a tau decay
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromTau(bool prompt_taus_only=false) const;

    /// @brief Determine whether the particle is from a prompt tau decay
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromPromptTau() const { return fromTau(true); }

    /// @brief Determine whether the particle is from a hadron or tau decay
    ///
    /// Specifically, walk up the ancestor chain until a status 2 hadron or
    /// tau is found, if at all.
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromDecay() const { return fromHadron() || fromPromptTau(); }

    /// @brief Shorthand definition of 'promptness' based on set definition flags
    ///
    /// @note This one doesn't make any judgements about final-stateness
    bool isPrompt(bool from_prompt_tau=false, bool from_prompt_mu=false) const;

    //@}


    /// @name Decay info
    //@{

    /// Whether this particle is stable according to the generator
    bool isStable() const;

    /// Get a list of the direct descendants from the current particle
    Particles children(const Cut& c=Cuts::open()) const;

    /// Get a list of all the descendants (including duplication of parents and children) from the current particle
    Particles allDescendants(const Cut& c=Cuts::open(), bool remove_duplicates=true) const;

    /// Get a list of all the stable descendants from the current particle
    /// @todo Use recursion through replica-avoiding MCUtils functions to avoid bookkeeping duplicates
    /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
    Particles stableDescendants(const Cut& c=Cuts::open()) const;

    /// Flight length (divide by mm or cm to get the appropriate units)
    double flightLength() const;

    //@}


  private:

    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;

    /// The PDG ID code for this Particle.
    PdgId _id;

    /// The momentum of this particle.
    FourMomentum _momentum;

    /// The creation position of this particle.
    FourVector _origin;

  };


  /// @name String representation and streaming support
  //@{

  /// Represent a Particle as a string.
  std::string to_str(const Particle& p);

  /// Allow a Particle to be passed to an ostream.
  inline std::ostream& operator<<(std::ostream& os, const Particle& p) {
    os << to_str(p);
    return os;
  }


  /// Represent a ParticlePair as a string.
  std::string to_str(const ParticlePair& pair);

  /// Allow ParticlePair to be passed to an ostream.
  inline std::ostream& operator<<(std::ostream& os, const ParticlePair& pp) {
    os << to_str(pp);
    return os;
  }

  //@}


}


#include "Rivet/Tools/ParticleUtils.hh"

#endif
