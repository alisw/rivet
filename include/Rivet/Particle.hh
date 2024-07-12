// -*- C++ -*-
#ifndef RIVET_Particle_HH
#define RIVET_Particle_HH

#include "Rivet/Particle.fhh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/RivetFastJet.hh"
#include "Rivet/Math/LorentzTrans.hh"
// NOTE: Rivet/Tools/ParticleUtils.hh included at the end

namespace Rivet {


  /// @brief Specialised vector of Particle objects.
  ///
  /// A specialised version of vector<Particle> which is able to implicitly and
  /// explicitly convert to a vector of FourMomentum.
  ///
  /// @todo Add explicit and implicit conversion to PseudoJets
  ///
  // typedef std::vector<Particle> Particles;
  class Particles : public std::vector<Particle> {
  public:
    using base = std::vector<Particle>; //< using-declarations don't like template syntax
    using base::base; //< import base-class constructors
    Particles();
    Particles(const std::vector<Particle>& vps);
    FourMomenta moms() const;
    PseudoJets pseudojets() const;
    operator FourMomenta () const { return moms(); }
    operator PseudoJets () const { return pseudojets(); }
    Particles& operator += (const Particle& p);
    Particles& operator += (const Particles& ps);
  };

  Particles operator + (const Particles& a, const Particles& b);

  /// Typedef for a pair of Particle objects.
  typedef std::pair<Particle, Particle> ParticlePair;

  //@}



  /////////////////////



  /// Particle representation, either from a HepMC::GenEvent or reconstructed.
  class Particle : public ParticleBase {
  public:

    /// @name Constructors
    /// @{

    /// Default constructor.
    /// @note A particle without info is useless. This only exists to keep STL containers happy.
    Particle()
      : ParticleBase(),
        _original(nullptr), _id(PID::ANY), _isDirect(4, std::make_pair(false,false))
    {   }

    /// Constructor from PID and momentum.
    Particle(PdgId pid, const FourMomentum& mom, const FourVector& pos=FourVector(), ConstGenParticlePtr gp=nullptr)
      : ParticleBase(),
        _original(gp), _id(pid),
        _momentum(mom), _origin(pos),
        _isDirect(4, std::make_pair(false,false))
    {   }

    /// Constructor from PID, momentum, and a GenParticle for relational links.
    Particle(PdgId pid, const FourMomentum& mom, ConstGenParticlePtr gp, const FourVector& pos=FourVector())
      : Particle(pid, mom, pos, gp)
    {   }

    /// Constructor from a HepMC GenParticle pointer.
    Particle(ConstGenParticlePtr gp)
      : ParticleBase(),
        _original(gp), _id(gp->pdg_id()),
        _momentum(gp->momentum()),
        _isDirect(4, std::make_pair(false,false))
    {
      ConstGenVertexPtr vprod = gp->production_vertex();
      if (vprod != nullptr) {
        setOrigin(vprod->position().t(), vprod->position().x(), vprod->position().y(), vprod->position().z());
      }
    }

    /// Constructor from a HepMC GenParticle reference.
    Particle(const RivetHepMC::GenParticle& gp)
      : Particle(HepMCUtils::getParticlePtr(gp))
    {   }

    /// @}

    /// @name Kinematic properties
    /// @{

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
    /// @{

    /// The origin position (and time).
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

    /// @}

    /// @name Displacement-projection properties
    /// @{

    /// Find the point of closest approach to the primary vertex
    Vector3 closestApproach() const {
      const FourVector& v0 = origin();
      /// @todo Check that this works with all angles
      const double rho0 = origin().perp() / sin(this->phi() - origin().phi());
      const double phi0 = M_PI/2 - this->phi();
      const double x0 = rho0 * cos(phi0);
      const double y0 = rho0 * sin(phi0);
      const double z0 = origin().z() - v0.perp()/tan(this->theta());
      return Vector3(x0, y0, z0);
    }

    /// @}


    /// @name Other representations and implicit casts to momentum-like objects
    /// @{

    /// Converter to FastJet3 PseudoJet
    virtual fastjet::PseudoJet pseudojet() const {
      return fastjet::PseudoJet(mom().px(), mom().py(), mom().pz(), mom().E());
    }

    /// Cast operator to FastJet3 PseudoJet
    operator PseudoJet () const { return pseudojet(); }


    /// Set a const pointer to the original GenParticle
    Particle& setGenParticle(ConstGenParticlePtr gp) {
      _original = gp;
      return *this;
    }

    /// Get a const pointer to the original GenParticle
    ConstGenParticlePtr genParticle() const {
      return _original;
    }

    /// Cast operator for conversion to GenParticle*
    /// @note Not implicit since that would enable accidental Particle::operator== comparisons
    explicit operator ConstGenParticlePtr () const { return genParticle(); }

    /// @}


    /// @name Particle ID code accessors
    /// @{

    /// This Particle's PDG ID code.
    PdgId pid() const { return _id; }
    /// Absolute value of the PDG ID code.
    PdgId abspid() const { return std::abs(_id); }

    /// @}


    /// @name Charge
    /// @{

    /// The charge of this Particle.
    double charge() const { return PID::charge(pid()); }

    /// The absolute charge of this Particle.
    double abscharge() const { return PID::abscharge(pid()); }

    /// Three times the charge of this Particle (i.e. integer multiple of smallest quark charge).
    int charge3() const { return PID::charge3(pid()); }

    /// Three times the absolute charge of this Particle (i.e. integer multiple of smallest quark charge).
    int abscharge3() const { return PID::abscharge3(pid()); }

    /// Is this Particle charged?
    bool isCharged() const { return charge3() != 0; }

    /// @}


    /// @name Particle species
    /// @{

    /// Is this a hadron?
    bool isHadron() const { return PID::isHadron(pid()); }

    /// Is this a meson?
    bool isMeson() const { return PID::isMeson(pid()); }

    /// Is this a baryon?
    bool isBaryon() const { return PID::isBaryon(pid()); }

    /// Is this a lepton?
    bool isLepton() const { return PID::isLepton(pid()); }

    /// Is this a charged lepton?
    bool isChargedLepton() const { return PID::isChargedLepton(pid()); }

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

    /// Is this a parton? (Hopefully not very often... fiducial FTW)
    bool isParton() const { return PID::isParton(pid()); }

    /// @}


    /// @name Constituents (for composite particles)
    /// @{

    /// Set direct constituents of this particle
    virtual void setConstituents(const Particles& cs, bool setmom=false);

    /// Add a single direct constituent to this particle
    virtual void addConstituent(const Particle& c, bool addmom=false);

    /// Add direct constituents to this particle
    virtual void addConstituents(const Particles& cs, bool addmom=false);


    /// Determine if this Particle is a composite of other Rivet Particles
    bool isComposite() const { return !constituents().empty(); }


    /// @brief Direct constituents of this particle, returned by reference
    ///
    /// The returned vector will be empty if this particle is non-composite,
    /// and its entries may themselves be composites.
    const Particles& constituents() const { return _constituents; }

    /// @brief Direct constituents of this particle, sorted by a functor
    /// @note Returns a copy, thanks to the sorting
    const Particles constituents(const ParticleSorter& sorter) const {
      return sortBy(constituents(), sorter);
    }

    /// @brief Direct constituents of this particle, filtered by a Cut
    /// @note Returns a copy, thanks to the filtering
    const Particles constituents(const Cut& c) const {
      return filter_select(constituents(), c);
    }

    /// @brief Direct constituents of this particle, sorted by a functor
    /// @note Returns a copy, thanks to the filtering and sorting
    const Particles constituents(const Cut& c, const ParticleSorter& sorter) const {
      return sortBy(constituents(c), sorter);
    }

    /// @brief Direct constituents of this particle, filtered by a selection functor
    /// @note Returns a copy, thanks to the filtering
    const Particles constituents(const ParticleSelector& selector) const {
      return filter_select(constituents(), selector);
    }

    /// @brief Direct constituents of this particle, filtered and sorted by functors
    /// @note Returns a copy, thanks to the filtering and sorting
    const Particles constituents(const ParticleSelector& selector, const ParticleSorter& sorter) const {
      return sortBy(constituents(selector), sorter);
    }


    /// @brief Fundamental constituents of this particle
    /// @note Returns {{*this}} if this particle is non-composite.
    Particles rawConstituents() const;

    /// @brief Fundamental constituents of this particle, sorted by a functor
    /// @note Returns a copy, thanks to the sorting
    const Particles rawConstituents(const ParticleSorter& sorter) const {
      return sortBy(rawConstituents(), sorter);
    }

    /// @brief Fundamental constituents of this particle, filtered by a Cut
    /// @note Returns a copy, thanks to the filtering
    const Particles rawConstituents(const Cut& c) const {
      return filter_select(rawConstituents(), c);
    }

    /// @brief Fundamental constituents of this particle, sorted by a functor
    /// @note Returns a copy, thanks to the filtering and sorting
    const Particles rawConstituents(const Cut& c, const ParticleSorter& sorter) const {
      return sortBy(rawConstituents(c), sorter);
    }

    /// @brief Fundamental constituents of this particle, filtered by a selection functor
    /// @note Returns a copy, thanks to the filtering
    const Particles rawConstituents(const ParticleSelector& selector) const {
      return filter_select(rawConstituents(), selector);
    }

    /// @brief Fundamental constituents of this particle, filtered and sorted by functors
    /// @note Returns a copy, thanks to the filtering and sorting
    const Particles rawConstituents(const ParticleSelector& selector, const ParticleSorter& sorter) const {
      return sortBy(rawConstituents(selector), sorter);
    }

    /// @}


    /// @name Ancestry (for fundamental particles with a HepMC link)
    /// @{

    /// Get a list of the direct parents of the current particle (with optional selection Cut)
    ///
    /// @note This is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    Particles parents(const Cut& c=Cuts::OPEN) const;

    /// Get a list of the direct parents of the current particle (with selector function)
    ///
    /// @note This is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    Particles parents(const ParticleSelector& f) const {
      return filter_select(parents(), f);
    }

    /// Check whether any particle in the particle's parent list has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasParentWith(const ParticleSelector& f) const {
      return !parents(f).empty();
    }
    /// Check whether any particle in the particle's parent list has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasParentWith(const Cut& c) const;

    /// Check whether any particle in the particle's parent list does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasParentWithout(const ParticleSelector& f) const {
      return hasParentWith([&](const Particle& p){ return !f(p); });
    }
    /// Check whether any particle in the particle's parent list does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasParentWithout(const Cut& c) const;

    /// Check whether a given PID is found in the particle's parent list
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    ///
    /// @deprecated Prefer e.g. hasParentWith(Cut::pid == 123)
    bool hasParent(PdgId pid) const;



    /// Get a list of the ancestors of the current particle (with optional selection Cut)
    ///
    /// @note By default only physical ancestors, with status=2, are returned.
    ///
    /// @note This is valid in MC, but may not be answerable experimentally --
    /// use this function with care when replicating experimental analyses!
    Particles ancestors(const Cut& c=Cuts::OPEN, bool only_physical=true) const;

    /// Get a list of the direct parents of the current particle (with selector function)
    ///
    /// @note By default only physical ancestors, with status=2, are returned.
    ///
    /// @note This is valid in MC, but may not be answerable experimentally --
    /// use this function with care when replicating experimental analyses!
    Particles ancestors(const ParticleSelector& f, bool only_physical=true) const {
      return filter_select(ancestors(Cuts::OPEN, only_physical), f);
    }

    /// Check whether any particle in the particle's ancestor list has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasAncestorWith(const ParticleSelector& f, bool only_physical=true) const {
      return !ancestors(f, only_physical).empty();
    }
    /// Check whether any particle in the particle's ancestor list has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasAncestorWith(const Cut& c, bool only_physical=true) const;

    /// Check whether any particle in the particle's ancestor list does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasAncestorWithout(const ParticleSelector& f, bool only_physical=true) const {
      return hasAncestorWith([&](const Particle& p){ return !f(p); }, only_physical);
    }
    /// Check whether any particle in the particle's ancestor list does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasAncestorWithout(const Cut& c, bool only_physical=true) const;

    /// Check whether a given PID is found in the particle's ancestor list
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    ///
    /// @deprecated Prefer hasAncestorWith(Cuts::pid == pid) etc.
    bool hasAncestor(PdgId pid, bool only_physical=true) const;


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

    /// @brief Determine whether the particle is from a tau which decayed hadronically
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool fromHadronicTau(bool prompt_taus_only=false) const;

    /// @brief Determine whether the particle is from a hadron or tau decay
    ///
    /// Specifically, walk up the ancestor chain until a status 2 hadron or
    /// tau is found, if at all.
    ///
    /// @note This question is valid in MC, but may not be perfectly answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    ///
    DEPRECATED("Too vague: use fromHadron() || fromPromptTau(), or isDirect()")
    bool fromDecay() const { return fromHadron() || fromPromptTau(); }

    /// @brief Shorthand definition of 'promptness' based on set definition flags
    ///
    /// A "direct" particle is one directly connected to the hard process. It is a
    /// preferred alias for "prompt", since it has no confusing implications about
    /// distinguishability by timing information.
    ///
    /// The boolean arguments allow a decay lepton to be considered direct if
    /// its parent was a "real" direct lepton.
    ///
    /// @note This one doesn't make any judgements about final-stateness
    bool isDirect(bool allow_from_direct_tau=false, bool allow_from_direct_mu=false) const;

    /// Alias for isDirect
    bool isPrompt(bool allow_from_prompt_tau=false, bool allow_from_prompt_mu=false) const {
      return isDirect(allow_from_prompt_tau, allow_from_prompt_mu);
    }

    /// @}


    /// @name Decay info
    /// @{

    /// Whether this particle is stable according to the generator
    bool isStable() const;

    /// @todo isDecayed? How to restrict to physical particles?


    /// Get a list of the direct descendants from the current particle (with optional selection Cut)
    Particles children(const Cut& c=Cuts::OPEN) const;

    /// Get a list of the direct descendants from the current particle (with selector function)
    Particles children(const ParticleSelector& f) const {
      return filter_select(children(), f);
    }

    /// Check whether any direct child of this particle has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasChildWith(const ParticleSelector& f) const {
      return !children(f).empty();
    }
    /// Check whether any direct child of this particle has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasChildWith(const Cut& c) const;

    /// Check whether any direct child of this particle does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasChildWithout(const ParticleSelector& f) const {
      return hasChildWith([&](const Particle& p){ return !f(p); });
    }
    /// Check whether any direct child of this particle does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasChildWithout(const Cut& c) const;


    /// Get a list of all the descendants from the current particle (with optional selection Cut)
    Particles allDescendants(const Cut& c=Cuts::OPEN, bool remove_duplicates=true) const;

    /// Get a list of all the descendants from the current particle (with selector function)
    Particles allDescendants(const ParticleSelector& f, bool remove_duplicates=true) const {
      return filter_select(allDescendants(Cuts::OPEN, remove_duplicates), f);
    }

    /// Check whether any descendant of this particle has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasDescendantWith(const ParticleSelector& f, bool remove_duplicates=true) const {
      return !allDescendants(f, remove_duplicates).empty();
    }
    /// Check whether any descendant of this particle has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasDescendantWith(const Cut& c, bool remove_duplicates=true) const;

    /// Check whether any descendant of this particle does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasDescendantWithout(const ParticleSelector& f, bool remove_duplicates=true) const {
      return hasDescendantWith([&](const Particle& p){ return !f(p); }, remove_duplicates);
    }
    /// Check whether any descendant of this particle does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasDescendantWithout(const Cut& c, bool remove_duplicates=true) const;


    /// Get a list of all the stable descendants from the current particle (with optional selection Cut)
    ///
    /// @todo Use recursion through replica-avoiding MCUtils functions to avoid bookkeeping duplicates
    /// @todo Insist that the current particle is post-hadronization, otherwise throw an exception?
    Particles stableDescendants(const Cut& c=Cuts::OPEN) const;

    /// Get a list of all the stable descendants from the current particle (with selector function)
    Particles stableDescendants(const ParticleSelector& f) const {
      return filter_select(stableDescendants(), f);
    }

    /// Check whether any stable descendant of this particle has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasStableDescendantWith(const ParticleSelector& f) const {
      return !stableDescendants(f).empty();
    }
    /// Check whether any stable descendant of this particle has the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasStableDescendantWith(const Cut& c) const;

    /// Check whether any stable descendant of this particle does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasStableDescendantWithout(const ParticleSelector& f) const {
      return hasStableDescendantWith([&](const Particle& p){ return !f(p); });
    }
    /// Check whether any stable descendant of this particle does not have the requested property
    ///
    /// @note This question is valid in MC, but may not be answerable
    /// experimentally -- use this function with care when replicating
    /// experimental analyses!
    bool hasStableDescendantWithout(const Cut& c) const;


    /// Flight length of the particle from origin to decay
    ///
    /// @note Divide by mm or cm as usual to get the appropriate units.
    double flightLength() const;

    /// @}


    /// @name Duplicate testing
    /// @{

    /// @brief Determine whether a particle is the first in a decay chain to meet the function requirement
    inline bool isFirstWith(const ParticleSelector& f) const {
      if (!f(*this)) return false; //< This doesn't even meet f, let alone being the last to do so
      if (any(parents(), f)) return false; //< If a direct parent has this property, this isn't the first
      return true;
    }

    /// @brief Determine whether a particle is the first in a decay chain not to meet the function requirement
    inline bool isFirstWithout(const ParticleSelector& f) const {
      return isFirstWith([&](const Particle& p){ return !f(p); });
    }

    /// @brief Determine whether a particle is the last in a decay chain to meet the function requirement
    inline bool isLastWith(const ParticleSelector& f) const {
      if (!f(*this)) return false; //< This doesn't even meet f, let alone being the last to do so
      if (any(children(), f)) return false; //< If a child has this property, this isn't the last
      return true;
    }

    /// @brief Determine whether a particle is the last in a decay chain not to meet the function requirement
    inline bool isLastWithout(const ParticleSelector& f) const {
      return isLastWith([&](const Particle& p){ return !f(p); });
    }

    /// @}


    /// @name Comparison
    /// @{

    /// Compare particles, based on "external" characteristics, with a little angular tolerance
    ///
    /// @note Not a deep comparison: GenParticle ptr and constituents are not used in the comparison
    bool isSame(const Particle& other) const {
      if (pid() != other.pid()) return false;
      if (!isZero((mom() - other.mom()).mod())) return false;
      if (!isZero((origin() - other.origin()).mod())) return false;
      return true;
    }

    /// @}


  protected:

    /// A pointer to the original GenParticle from which this Particle is projected.
    ConstGenParticlePtr _original;

    /// Constituent particles if this is a composite (may be empty)
    Particles _constituents;

    /// The PDG ID code for this Particle.
    PdgId _id;

    /// The momentum of this particle.
    FourMomentum _momentum;

    /// The creation position of this particle.
    FourVector _origin;

    /// Cached computation of directness, via ancestry. Second element is cache status
    /// @todo Replace this awkward caching with C++17 std::optional
    mutable std::vector<std::pair<bool,bool> > _isDirect;

  };


  /// @name String representation and streaming support
  /// @{

  /// Allow a Particle to be passed to an ostream.
  std::ostream& operator << (std::ostream& os, const Particle& p);

  /// Allow ParticlePair to be passed to an ostream.
  std::ostream& operator << (std::ostream& os, const ParticlePair& pp);

  /// @}


}


#include "Rivet/Tools/ParticleUtils.hh"

#endif
