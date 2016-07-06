// -*- C++ -*-
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  PromptFinalState::PromptFinalState(const FinalState& fsp)
    : _acceptMuDecays(false), _acceptTauDecays(false)
  {
    setName("PromptFinalState");
    addProjection(fsp, "FS");
  }


  PromptFinalState::PromptFinalState(const Cut & c)
    : _acceptMuDecays(false), _acceptTauDecays(false)
  {
    setName("PromptFinalState");
    addProjection(FinalState(c), "FS");
  }

  int PromptFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;
    const PromptFinalState& other = dynamic_cast<const PromptFinalState&>(p);
    return cmp(_acceptMuDecays, other._acceptMuDecays) || cmp(_acceptTauDecays, other._acceptTauDecays);
  }


  bool PromptFinalState::isPrompt(const Particle& p) const {
    if (p.genParticle() == NULL) return false; // no HepMC connection, give up! Throw UserError exception?
    /// @todo Shouldn't a const vertex be being returned? Ah, HepMC...
    GenVertex* prodVtx = p.genParticle()->production_vertex();
    if (prodVtx == NULL) return false; // orphaned particle, has to be assume false
    const pair<GenParticle*, GenParticle*> beams = prodVtx->parent_event()->beam_particles();

    /// @todo Would be nicer to be able to write this recursively up the chain, exiting as soon as a parton or string/cluster is seen
    foreach (const GenParticle* ancestor, Rivet::particles(prodVtx, HepMC::ancestors)) {
      const PdgId pid = ancestor->pdg_id();
      if (ancestor->status() != 2) continue; // no non-standard statuses or beams to be used in decision making
      if (ancestor == beams.first || ancestor == beams.second) continue; // PYTHIA6 uses status 2 for beams, I think... (sigh)
      if (PID::isParton(pid)) continue; // PYTHIA6 also uses status 2 for some partons, I think... (sigh)
      if (PID::isHadron(pid)) return false; // prompt particles can't be from hadron decays
      if (abs(pid) == PID::TAU && p.abspid() != PID::TAU && !_acceptTauDecays) return false; // allow or ban particles from tau decays (permitting tau copies)
      if (abs(pid) == PID::MUON && p.abspid() != PID::MUON && !_acceptMuDecays) return false; // allow or ban particles from muon decays (permitting muon copies)
    }
    return true;
  }


  void PromptFinalState::project(const Event& e) {
    _theParticles.clear();

    const Particles& particles = applyProjection<FinalState>(e, "FS").particles();
    foreach (const Particle& p, particles)
      if (isPrompt(p)) _theParticles.push_back(p);
    MSG_DEBUG("Number of final state particles not from hadron decays = " << _theParticles.size());

    if (getLog().isActive(Log::TRACE)) {
      foreach (const Particle& p, _theParticles)
        MSG_TRACE("Selected: " << p.pid() << ", charge = " << p.charge());
    }
  }


}
