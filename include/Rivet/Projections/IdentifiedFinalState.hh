// -*- C++ -*-
#ifndef RIVET_IdentifiedFinalState_HH
#define RIVET_IdentifiedFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Produce a final state which only contains specified particle IDs.
  class IdentifiedFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor with a FinalState and optional list of PDG ID codes.
    IdentifiedFinalState(const FinalState& fsp, const vector<PdgId>& pids=vector<PdgId>());

    /// Constructor with a list of PDG ID codes and a FinalState.
    IdentifiedFinalState(const vector<PdgId>& pids, const FinalState& fsp);

    /// Constructor with a FinalState and a single of PDG ID code.
    IdentifiedFinalState(const FinalState& fsp, PdgId pid);

    /// Constructor with a single PDG ID code and a FinalState.
    IdentifiedFinalState(PdgId pid, const FinalState& fsp);


    /// Construction using optional Cuts object and optional list of PDG ID codes
    IdentifiedFinalState(const Cut& c=Cuts::open(), const vector<PdgId>& pids=vector<PdgId>());

    /// Construction using list of PDG ID codes and an optional Cuts object
    IdentifiedFinalState(const vector<PdgId>& pids, const Cut& c=Cuts::open());

    /// Construction using Cuts object and a single PDG ID code
    IdentifiedFinalState(const Cut& c, PdgId pid);

    /// Construction using a single PDG ID code and an optional Cuts object
    IdentifiedFinalState(PdgId pid, const Cut& c=Cuts::open());


    /// Constructor with eta range and pT_min arguments and optional list of PDG ID codes.
    /// @deprecated Use the versions with Cut or FinalState arguments
    DEPRECATED("Use the versions with Cut or FinalState arguments.")
    IdentifiedFinalState(double etamin, double etamax, double ptMin=0.0*GeV);


    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new IdentifiedFinalState(*this);
    }

    //@}


  public:

    /// Get the list of particle IDs to accept.
    const set<PdgId>& acceptedIds() const {
      return _pids;
    }

    /// Add an accepted particle ID.
    IdentifiedFinalState& acceptId(PdgId pid) {
      _pids.insert(pid);
      return *this;
    }

    /// Add a set of accepted particle IDs.
    IdentifiedFinalState& acceptIds(const vector<PdgId>& pids) {
      foreach (const PdgId pid, pids) {
        _pids.insert(pid);
      }
      return *this;
    }

    /// Add an accepted particle ID and its antiparticle.
    IdentifiedFinalState& acceptIdPair(PdgId pid) {
      _pids.insert(pid);
      _pids.insert(-pid);
      return *this;
    }

    /// Add a set of accepted particle IDs and their antiparticles.
    IdentifiedFinalState& acceptIdPairs(const vector<PdgId>& pids) {
      foreach (const PdgId pid, pids) {
        _pids.insert(pid);
        _pids.insert(-pid);
      }
      return *this;
    }

    /// Accept all neutrinos (convenience method).
    IdentifiedFinalState& acceptNeutrinos() {
      acceptIdPair(PID::NU_E);
      acceptIdPair(PID::NU_MU);
      acceptIdPair(PID::NU_TAU);
      return *this;
    }

    /// Accept all charged leptons (convenience method).
    IdentifiedFinalState& acceptChLeptons() {
      acceptIdPair(PID::ELECTRON);
      acceptIdPair(PID::MUON);
      acceptIdPair(PID::TAU);
      return *this;
    }

    /// Reset the list of particle IDs to accept.
    void reset() {
      _pids.clear();
    }

    // The remaining particles
    virtual const Particles& remainingParticles() const {
      return _remainingParticles;
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    /// The final-state particles.
    set<PdgId> _pids;

    // A vector of all other particles in the final state
    Particles _remainingParticles;

  };


}


#endif
