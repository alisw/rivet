// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Particle.fhh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Utils.hh"
// NOTE: Rivet/Tools/ParticleUtils.hh included at the end
#include "fastjet/PseudoJet.hh"

namespace Rivet {


  /// Representation of particles from a HepMC::GenEvent.
  class Particle : public ParticleBase {
  public:

    /// Default constructor.
    /// @note A particle without info is useless. This only exists to keep STL containers happy.
    Particle()
      : ParticleBase(),
        _original(0), _id(0), _momentum()
    { }

    /// Constructor without GenParticle.
    Particle(PdgId pid, const FourMomentum& mom)
      : ParticleBase(),
        _original(0), _id(pid), _momentum(mom)
    { }

    /// Constructor from a HepMC GenParticle.
    Particle(const GenParticle& gp)
      : ParticleBase(),
        _original(&gp), _id(gp.pdg_id()),
        _momentum(gp.momentum())
    { }

    /// Constructor from a HepMC GenParticle pointer.
    Particle(const GenParticle* gp)
      : ParticleBase(),
        _original(gp), _id(gp->pdg_id()),
        _momentum(gp->momentum())
    { }


  public:

    /// Converter to FastJet3 PseudoJet
    virtual fastjet::PseudoJet pseudojet() const {
      return fastjet::PseudoJet(mom().px(), mom().py(), mom().pz(), mom().E());
    }

    /// Cast operator to FastJet3 PseudoJet
    operator PseudoJet () const { return pseudojet(); }


  public:

    /// @name Basic particle specific properties
    //@{

    /// Get a const reference to the original GenParticle.
    const GenParticle* genParticle() const {
      return _original;
    }

    /// Cast operator for conversion to GenParticle*
    operator const GenParticle* () const { return genParticle(); }

    /// The momentum.
    const FourMomentum& momentum() const {
      return _momentum;
    }
    /// Set the momentum.
    Particle& setMomentum(const FourMomentum& momentum) {
      _momentum = momentum;
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
    double charge() const {
      return PID::charge(pid());
    }
    /// Three times the charge of this Particle (i.e. integer multiple of smallest quark charge).
    int threeCharge() const {
      return PID::threeCharge(pid());
    }

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

    //@}


    /// @name Ancestry properties
    //@{

    /// Check whether a given PID is found in the GenParticle's ancestor list
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

    //@}


    /// @name Decay info
    //@{

    /// Whether this particle is stable according to the generator
    bool isStable() const {
      return genParticle() != NULL && genParticle()->status() == 1 && genParticle()->end_vertex() == NULL;
    }

    /// Get a list of the direct descendants from the current particle
    vector<Particle> children() const {
      vector<Particle> rtn;
      if (isStable()) return rtn;
      /// @todo Remove this const mess crap when HepMC doesn't suck
      HepMC::GenVertex* gv = const_cast<HepMC::GenVertex*>( genParticle()->end_vertex() );
      /// @todo Would like to do this, but the range objects are broken
      // foreach (const GenParticle* gp, gv->particles(HepMC::children))
      //   rtn += Particle(gp);
      for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::children); it != gv->particles_end(HepMC::children); ++it)
        rtn += Particle(*it);
      return rtn;
    }

    /// Get a list of all the descendants (including duplication of parents and children) from the current particle
    /// @todo Use recursion through replica-avoiding MCUtils functions to avoid bookkeeping duplicates
    /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
    vector<Particle> allDescendants() const {
      vector<Particle> rtn;
      if (isStable()) return rtn;
      /// @todo Remove this const mess crap when HepMC doesn't suck
      HepMC::GenVertex* gv = const_cast<HepMC::GenVertex*>( genParticle()->end_vertex() );
      /// @todo Would like to do this, but the range objects are broken
      // foreach (const GenParticle* gp, gv->particles(HepMC::descendants))
      //   rtn += Particle(gp);
      for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it != gv->particles_end(HepMC::descendants); ++it)
        rtn += Particle(*it);
      return rtn;
    }

    /// Get a list of all the stable descendants from the current particle
    /// @todo Use recursion through replica-avoiding MCUtils functions to avoid bookkeeping duplicates
    /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
    vector<Particle> stableDescendants() const {
      vector<Particle> rtn;
      if (isStable()) return rtn;
      /// @todo Remove this const mess crap when HepMC doesn't suck
      HepMC::GenVertex* gv = const_cast<HepMC::GenVertex*>( genParticle()->end_vertex() );
      /// @todo Would like to do this, but the range objects are broken
      // foreach (const GenParticle* gp, gv->particles(HepMC::descendants))
      //   if (gp->status() == 1 && gp->end_vertex() == NULL)
      //     rtn += Particle(gp);
      for (GenVertex::particle_iterator it = gv->particles_begin(HepMC::descendants); it != gv->particles_end(HepMC::descendants); ++it)
        if ((*it)->status() == 1 && (*it)->end_vertex() == NULL)
          rtn += Particle(*it);
      return rtn;
    }

    /// Flight length (divide by mm or cm to get the appropriate units)
    double flightLength() const {
      if (isStable()) return -1;
      if (genParticle() == NULL) return 0;
      if (genParticle()->production_vertex() == NULL) return 0;
      const HepMC::FourVector v1 = genParticle()->production_vertex()->position();
      const HepMC::FourVector v2 = genParticle()->end_vertex()->position();
      return sqrt(sqr(v2.x()-v1.x()) + sqr(v2.y()-v1.y()) + sqr(v2.z()-v1.z()));
    }

    //@}


  private:

    /// A pointer to the original GenParticle from which this Particle is projected.
    const GenParticle* _original;

    /// The PDG ID code for this Particle.
    PdgId _id;

    /// The momentum of this projection of the Particle.
    FourMomentum _momentum;

  };


  /// @name String representation
  //@{

  /// Print a ParticlePair as a string.
  inline std::string to_str(const ParticlePair& pair) {
    stringstream out;
    out << "["
        << PID::toParticleName(pair.first.pid()) << " @ "
        << pair.first.momentum().E()/GeV << " GeV, "
        << PID::toParticleName(pair.second.pid()) << " @ "
        << pair.second.momentum().E()/GeV << " GeV]";
    return out.str();
  }

  /// Allow ParticlePair to be passed to an ostream.
  inline std::ostream& operator<<(std::ostream& os, const ParticlePair& pp) {
    os << to_str(pp);
    return os;
  }

  //@}


}


#include "Rivet/Tools/ParticleUtils.hh"

#endif
