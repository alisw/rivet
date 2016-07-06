// -*- C++ -*-
#ifndef RIVET_PromptFinalState_HH
#define RIVET_PromptFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Find final state particles directly connected to the hard process.
  ///
  /// The definition of "prompt" used in Rivet is that from high-scale physics, i.e.
  /// particles directly connected to the hard process in an interaction, regardless
  /// of realistic reconstructibility of displaced vertices, etc. By construction
  /// hadrons cannot be considered prompt as they will be colour connected to other
  /// parts of the event through non-perturbative effects: this projection can
  /// return electrons, muons, photons, and exotic particles which do not have a
  /// hadron in their post-hadronization ancestor chain. Flags exist to choose
  /// whether intermediate tau or muon decays invalidate a particle's promptness.
  ///
  /// @todo Decide how to treat brem photons off prompt leptons -- are they also prompt? "Decay" does not change the lepton PID...
  class PromptFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{
    // Final State
    PromptFinalState(const FinalState& fsp);

    /// Cut constructor.
    PromptFinalState(const Cut & c);

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new PromptFinalState(*this);
    }
    //@}

    /// Accept particles from decays of prompt muons as themselves being prompt?
    void acceptMuonDecays(bool acc=true) { _acceptMuDecays = acc; }
    /// Accept particles from decays of prompt taus as themselves being prompt?
    void acceptTauDecays(bool acc=true) { _acceptTauDecays = acc; }

    /// Decide if a given particle is prompt based on set definition flags
    /// @todo Move into ParticleUtils / MCUtils
    /// @note This one doesn't make any judgements about final-stateness
    bool isPrompt(const Particle& p) const;

  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;

  private:

    bool _acceptMuDecays, _acceptTauDecays;

  };

}


#endif
