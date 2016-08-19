// -*- C++ -*-
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @todo Reduce the cut & paste duplication between the constructors. With C++11 constructors can chain...


  DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                                 double dRmax, const Cut& cut, bool cluster, bool useDecayPhotons)
    : FinalState(cut),
      _dRmax(dRmax), _cluster(cluster), _fromDecay(useDecayPhotons)
  {
    setName("DressedLeptons");
    IdentifiedFinalState photonfs(photons);
    photonfs.acceptId(PID::PHOTON);
    addProjection(photonfs, "Photons");
    addProjection(bareleptons, "Leptons");
  }


  DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                                 double dRmax, bool cluster, const Cut& cut,
                                 bool useDecayPhotons)
    : FinalState(cut),
      _dRmax(dRmax), _cluster(cluster), _fromDecay(useDecayPhotons)
  {
    setName("DressedLeptons");
    IdentifiedFinalState photonfs(photons);
    photonfs.acceptId(PID::PHOTON);
    addProjection(photonfs, "Photons");
    addProjection(bareleptons, "Leptons");
  }


  DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                                 double dRmax, bool cluster,
                                 double etaMin, double etaMax,
                                 double pTmin, bool useDecayPhotons)
    : FinalState(etaMin, etaMax, pTmin),
      _dRmax(dRmax), _cluster(cluster), _fromDecay(useDecayPhotons)
  {
    setName("DressedLeptons");
    IdentifiedFinalState photonfs(photons);
    photonfs.acceptId(PID::PHOTON);
    addProjection(photonfs, "Photons");
    addProjection(bareleptons, "Leptons");
  }



  int DressedLeptons::compare(const Projection& p) const {
    // Compare the two as final states (for pT and eta cuts)
    const DressedLeptons& other = dynamic_cast<const DressedLeptons&>(p);
    int fscmp = FinalState::compare(other);
    if (fscmp != EQUIVALENT) return fscmp;

    const PCmp phcmp = mkNamedPCmp(p, "Photons");
    if (phcmp != EQUIVALENT) return phcmp;

    const PCmp sigcmp = mkNamedPCmp(p, "Leptons");
    if (sigcmp != EQUIVALENT) return sigcmp;

    return (cmp(_dRmax, other._dRmax) ||
            cmp(_cluster, other._cluster) ||
            cmp(_fromDecay, other._fromDecay));
  }


  void DressedLeptons::project(const Event& e) {
    _theParticles.clear();
    _clusteredLeptons.clear();

    const FinalState& signal = applyProjection<FinalState>(e, "Leptons");
    Particles bareleptons = signal.particles();
    if (bareleptons.empty()) return;

    vector<DressedLepton> allClusteredLeptons;
    for (size_t i = 0; i < bareleptons.size(); ++i) {
      allClusteredLeptons.push_back(DressedLepton(bareleptons[i]));
    }

    // Match each photon to its closest charged lepton within the dR cone
    const FinalState& photons = applyProjection<FinalState>(e, "Photons");
    for (const Particle& photon : photons.particles()) {
      // Ignore photon if it's from a hadron/tau decay and we're avoiding those
      if (!_fromDecay && photon.fromDecay()) continue;
      const FourMomentum p_P = photon.momentum();
      double dRmin = _dRmax;
      int idx = -1;
      for (size_t i = 0; i < bareleptons.size(); ++i) {
        // Only cluster photons around *charged* signal particles
        if (PID::threeCharge(bareleptons[i].pid()) == 0) continue;
        // Find the closest lepton
        const FourMomentum& p_l = bareleptons[i].momentum();
        double dR = deltaR(p_l, p_P);
        if (dR < dRmin) {
          dRmin = dR;
          idx = i;
        }
      }
      if (idx > -1) {
        if (_cluster) allClusteredLeptons[idx].addPhoton(photon, _cluster);
      }
    }

    for (const DressedLepton& lepton : allClusteredLeptons) {
      if (accept(lepton)) {
        _clusteredLeptons.push_back(lepton);
        _theParticles.push_back(lepton.constituentLepton());
        /// @todo Can't we use += here?
        _theParticles.insert(_theParticles.end(),
                             lepton.constituentPhotons().begin(),
                             lepton.constituentPhotons().end());
      }
    }
  }
}
