// -*- C++ -*-
#ifndef RIVET_Jet_HH
#define RIVET_Jet_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Jet.fhh"
#include "Rivet/Particle.hh"
#include "Rivet/Cuts.hh"
#include "Rivet/Tools/ParticleUtils.hh"
#include "fastjet/PseudoJet.hh"
#include <numeric>

namespace Rivet {


  /// @brief Representation of a clustered jet of particles.
  class Jet : public ParticleBase {
  public:

    /// @name Constructors
    //@{

    /// Constructor from a FastJet PseudoJet, with optional full particle constituents information.
    Jet(const fastjet::PseudoJet& pj, const Particles& particles=Particles(), const Particles& tags=Particles()) {
      setState(pj, particles, tags);
    }

    /// Set the jet data, with optional full particle information.
    Jet(const FourMomentum& pjet, const Particles& particles=Particles(), const Particles& tags=Particles()) {
      setState(pjet, particles, tags);
    }

    /// Set all the jet data, with full particle information.
    /// @deprecated Prefer the form where the 4-vec comes first and the particles list is optional.
    DEPRECATED("Prefer the form where the 4-vec comes first and the particles list is optional.")
    Jet(const Particles& particles, const FourMomentum& pjet) {
      setState(pjet, particles);
    }

    /// Default constructor -- only for STL storability
    Jet() { clear(); }

    //@}


    /// @name Access jet constituents
    //@{

    /// Number of particles in this jet.
    size_t size() const { return _particles.size(); }

    /// Get the particles in this jet.
    vector<Particle>& particles() { return _particles; }
    /// Get the particles in this jet (const version)
    const vector<Particle>& particles() const { return _particles; }

    /// Get the particles in this jet (FastJet-like alias)
    vector<Particle>& constituents() { return particles(); }
    /// Get the particles in this jet (FastJet-like alias, const version)
    const vector<Particle>& constituents() const { return particles(); }

    /// Check whether this jet contains a particular particle.
    bool containsParticle(const Particle& particle) const;
    /// Nicer alias for containsParticleId
    bool containsPID(const Particle& particle) const { return containsParticle(particle); }

    /// Check whether this jet contains a certain particle type.
    bool containsParticleId(PdgId pid) const;
    /// Nicer alias for containsParticleId
    bool containsPID(PdgId pid) const { return containsParticleId(pid); }

    /// Check whether this jet contains at least one of certain particle types.
    bool containsParticleId(const vector<PdgId>& pids) const;
    /// Nicer alias for containsParticleId
    bool containsPID(const vector<PdgId>& pids) const { return containsParticleId(pids); }

    //@}


    /// @name Tagging
    ///
    /// @note General sources of tag particles are planned. The default jet finding
    /// adds b-hadron, c-hadron, and tau tags by ghost association.
    //@{

    /// @brief Particles which have been tag-matched to this jet
    Particles& tags() { return _tags; }
    /// @brief Particles which have been tag-matched to this jet (const version)
    const Particles& tags() const { return _tags; }
    /// @brief Particles which have been tag-matched to this jet _and_ pass a Cut
    ///
    /// @note Note the less efficient return by value, due to the cut-pass filtering.
    Particles tags(const Cut& c) const;


    /// @brief b particles which have been tag-matched to this jet (and pass an optional Cut)
    ///
    /// The default jet finding adds b-hadron tags by ghost association.
    Particles bTags(const Cut& c=Cuts::open()) const;
    /// Does this jet have at least one b-tag (that passes an optional Cut)?
    bool bTagged(const Cut& c=Cuts::open()) const { return !bTags(c).empty(); }


    /// @brief c (and not b) particles which have been tag-matched to this jet (and pass an optional Cut)
    ///
    /// The default jet finding adds c-hadron tags by ghost association.
    Particles cTags(const Cut& c=Cuts::open()) const;
    /// Does this jet have at least one c-tag (that passes an optional Cut)?
    bool cTagged(const Cut& c=Cuts::open()) const { return !cTags(c).empty(); }


    /// @brief Tau particles which have been tag-matched to this jet (and pass an optional Cut)
    ///
    /// The default jet finding adds tau tags by ghost association.
    Particles tauTags(const Cut& c=Cuts::open()) const;
    /// Does this jet have at least one tau-tag (that passes an optional Cut)?
    bool tauTagged(const Cut& c=Cuts::open()) const { return !tauTags(c).empty(); }


    /// @brief Check whether this jet contains a bottom-flavoured hadron.
    ///
    /// @deprecated The bTags() or bTagged() function is probably what you want
    /// for tagging. This one ignores the tags() list and draws conclusions
    /// based directly on the jet constituents; the other gives a much better match
    /// to typical experimental methods.
    ///
    /// @note The decision is made by first trying to find a bottom-flavoured particle
    /// in the particles list. Most likely this will fail unless bottom hadrons
    /// are set stable. If @a include_decay_products is true (the default), a
    /// fallback is attempted, using the post-hadronization ancestor history of
    /// all constituents.
    //DEPRECATED("Prefer the bTags() or bTagged() function")
    bool containsBottom(bool include_decay_products=true) const;

    /// @brief Check whether this jet contains a charm-flavoured hadron.
    ///
    /// @deprecated The cTags() or cTagged() function is probably what you want
    /// for tagging. This one ignores the tags() list and draws conclusions
    /// based directly on the jet constituents; the other gives a much better match
    /// to typical experimental methods.
    ///
    /// @note The decision is made by first trying to find a charm-flavoured particle
    /// in the particles list. Most likely this will fail unless charmed hadrons
    /// are set stable. If @a include_decay_products is true (the default), a
    /// fallback is attempted, using the post-hadronization ancestor history of
    /// all constituents.
    //DEPRECATED("Prefer the cTags() or cTagged() function")
    bool containsCharm(bool include_decay_products=true) const;

    //@}


    /// @name Access additional effective jet 4-vector properties
    //@{

    /// Get equivalent single momentum four-vector.
    const FourMomentum& momentum() const { return _momentum; }

    /// Get the total energy of this jet.
    double totalEnergy() const { return momentum().E(); }

    /// Get the energy carried in this jet by neutral particles.
    double neutralEnergy() const;

    /// Get the energy carried in this jet by hadrons.
    double hadronicEnergy() const;

    //@}


    /// @name Interaction with FastJet
    //@{

    /// Access the internal FastJet3 PseudoJet (as a const reference)
    const fastjet::PseudoJet& pseudojet() const { return _pseudojet; }

    /// Cast operator to FastJet3 PseudoJet (as a const reference)
    operator const fastjet::PseudoJet& () const { return pseudojet(); }

    //@}


    /// @name Set the jet constituents and properties
    //@{

    /// @brief Set the jet data from a FastJet PseudoJet, with optional particle constituents and tags lists.
    ///
    /// @note The particles() list will be extracted from PseudoJet constituents
    /// by default, making use of an attached user info if one is found.
    Jet& setState(const fastjet::PseudoJet& pj, const Particles& particles=Particles(), const Particles& tags=Particles());

    /// Set all the jet data, with optional full particle constituent and tag information.
    Jet& setState(const FourMomentum& mom, const Particles& particles, const Particles& tags=Particles());

    /// @deprecated Prefer the 4-mom first-arg versions
    DEPRECATED("Prefer the 4-mom first-arg versions")
    Jet& setState(const Particles& particles, const FourMomentum& mom) { return setState(mom, particles); }

    /// @brief Set the particles collection with full particle information.
    ///
    /// If set, this overrides particle info extracted from the PseudoJet
    Jet& setParticles(const Particles& particles);
    Jet& setConstituents(const Particles& particles) { return setParticles(particles); }

    /// Reset this jet as empty.
    Jet& clear();

    //@}


  private:

    /// FJ3 PseudoJet member to unify PseudoJet and Jet
    fastjet::PseudoJet _pseudojet;

    /// Full constituent particle information. (Filled from PseudoJet if possible.)
    /// @todo Make these mutable or similar? Add a flag to force a cache rebuild?
    Particles _particles;

    /// Particles used to tag this jet (can be anything, but c and b hadrons are the most common)
    Particles _tags;

    /// Effective jet 4-vector (just for caching)
    mutable FourMomentum _momentum;

  };


}

#endif
