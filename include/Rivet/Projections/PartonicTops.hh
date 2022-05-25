// -*- C++ -*-
#ifndef RIVET_PartonicTops_HH
#define RIVET_PartonicTops_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Convenience finder of partonic top quarks
  ///
  /// @warning This projection requires there to be tops in the event record:
  /// there is no guarantee of this, especially where the top quark is treated
  /// (correctly) as a resonance rather than on-shell. Further, there is no
  /// guarantee that the kinematics assigned to such tops are consistent,
  /// physical, or even associated with the lab frame. A fiducial pseudo-top
  /// analysis approach is *strongly* recommended instead.
  class PartonicTops : public ParticleFinder {
  public:


    /// @brief Enum for categorising top quark decay modes
    ///
    /// More specifically, the decay mode of the W from the top. We presume top decay to a W and b quark.
    enum class DecayMode {
      ANY = 0,
      ALL = 0,
      ELECTRON,
      MUON,
      TAU,
      E_MU,
      E_MU_TAU,
      HADRONIC
    };

    /// @brief Enum for categorising which top quark to be selected: last (weakly decaying) or first?
    enum class WhichTop { FIRST, LAST };


    /// @name Constructors
    //@{

    /// Constructor taking decay mode details (and an optional cuts object)
    PartonicTops(DecayMode decaymode, bool emu_from_prompt_tau=true, bool include_hadronic_taus=false, const Cut& c=Cuts::OPEN, WhichTop whichtop=WhichTop::LAST)
      : ParticleFinder(c), _topmode(whichtop), _decaymode(decaymode),
        _emu_from_prompt_tau(emu_from_prompt_tau), _include_hadronic_taus(include_hadronic_taus)
    {  }

    /// Constructor taking decay mode details (and a non-optional cuts object)
    PartonicTops(DecayMode decaymode, const Cut& c, bool emu_from_prompt_tau=true, bool include_hadronic_taus=false, WhichTop whichtop=WhichTop::LAST)
      : PartonicTops(decaymode, emu_from_prompt_tau, include_hadronic_taus, c, whichtop)
    {  }

    /// Simple constructor optionally taking cuts object
    PartonicTops(const Cut& c=Cuts::OPEN, WhichTop whichtop=WhichTop::LAST)
      : PartonicTops(DecayMode::ALL, true, false, c, whichtop)
    {  }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PartonicTops);

    //@}


    /// Access to the found partonic tops
    const Particles& tops() const { return _theParticles; }


    /// Clear the projection
    void clear() {
      _theParticles.clear();
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& event) {

      // Warn about how terrible this is, the first time it's called!
      static bool donerubric = false;
      if (!donerubric) {
        MSG_WARNING("PartonicTops is not recommended: MC generators do not guarantee physical properties for, or even the existence of, partonic event-record entries. Caveat emptor!");
        donerubric = true;
      }

      // Find partonic tops
      _theParticles = filter_select(event.allParticles(_cuts), (_topmode == WhichTop::LAST ? lastParticleWith(isTop) : firstParticleWith(isTop)));

      // Filtering by decay mode
      if (_decaymode != DecayMode::ALL) {
        const auto decaycheck = [&](const Particle& t) {
          const Particles descendants = t.allDescendants();
          const bool prompt_e = any(descendants, [&](const Particle& p){ return p.abspid() == PID::ELECTRON && p.isPrompt(_emu_from_prompt_tau) && !p.hasAncestor(PID::PHOTON, false); });
          const bool prompt_mu = any(descendants, [&](const Particle& p){ return p.abspid() == PID::MUON && p.isPrompt(_emu_from_prompt_tau) && !p.hasAncestor(PID::PHOTON, false); });
          if (prompt_e && (_decaymode == DecayMode::ELECTRON || _decaymode == DecayMode::E_MU || _decaymode == DecayMode::E_MU_TAU)) return true;
          if (prompt_mu && (_decaymode == DecayMode::MUON || _decaymode == DecayMode::E_MU || _decaymode == DecayMode::E_MU_TAU)) return true;
          const bool prompt_tau = any(descendants, [&](const Particle& p){ return p.abspid() == PID::TAU && p.isPrompt()  && !p.hasAncestor(PID::PHOTON, false); });
          const bool prompt_hadronic_tau = any(descendants, [&](const Particle& p){ return p.abspid() == PID::TAU && p.isPrompt() && !p.hasAncestor(PID::PHOTON, false) && none(p.children(), isChargedLepton); });
          if (prompt_tau && (_decaymode == DecayMode::TAU || _decaymode == DecayMode::E_MU_TAU)) return (_include_hadronic_taus || !prompt_hadronic_tau);
          if (_decaymode == DecayMode::HADRONIC && (!prompt_e && !prompt_mu && (!prompt_tau || (_include_hadronic_taus && prompt_hadronic_tau)))) return true; //< logical hairiness...
          return false;
        };
        ifilter_select(_theParticles, decaycheck);
      }

      // Filtering and warning about unphysical partonic tops
      const auto physcheck = [&](const Particle& t) {
        if (t.E() < 0 || t.mass() < 0) {
          MSG_WARNING("Unphysical partonic top with negative E or m found: " << t.mom());
          return false;
        }
        return true;
      };
      ifilter_select(_theParticles, physcheck);
    }


    /// Compare projections.
    CmpState compare(const Projection& p) const {
      const PartonicTops& other = dynamic_cast<const PartonicTops&>(p);
      return cmp(_cuts, other._cuts) ||
        cmp(_topmode, other._topmode) ||
        cmp(_decaymode, other._decaymode) ||
        cmp(_emu_from_prompt_tau, other._emu_from_prompt_tau) ||
        cmp(_include_hadronic_taus, other._include_hadronic_taus);
    }


  private:

    WhichTop _topmode;

    DecayMode _decaymode;

    bool _emu_from_prompt_tau, _include_hadronic_taus;

  };


}


#endif
