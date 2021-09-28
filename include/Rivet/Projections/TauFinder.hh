#ifndef RIVET_TauFinder_HH
#define RIVET_TauFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Convenience finder of unstable taus
  ///
  /// @todo Convert to a general ParticleFinder, since it's not a true final state? Needs some care...
  class TauFinder : public FinalState {
  public:

    enum class DecayMode {
      ANY = 0,
      ALL = 0,
      LEPTONIC,
      HADRONIC
    };

    static bool isHadronic(const Particle& tau) {
      assert(tau.abspid() == PID::TAU);
      return any(tau.stableDescendants(), isHadron);
    }

    static bool isLeptonic(const Particle& tau) {
      return !isHadronic(tau);
    }


    TauFinder(DecayMode decaymode=DecayMode::ANY, const Cut& cut=Cuts::open()) {
      /// @todo What about directness/promptness?
      setName("TauFinder");
      _decmode = decaymode;
      declare(UnstableParticles(cut), "UFS");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(TauFinder);


    const Particles& taus() const { return _theParticles; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare with other projections.
    virtual CmpState compare(const Projection& p) const;


  private:

    /// The decaymode enum
    DecayMode _decmode;

  };


  /// @todo Make this the canonical name in future
  using Taus = TauFinder;


}


#endif
