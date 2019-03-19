// -*- C++ -*-
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  // On DressedLepton helper class
  //{

  DressedLepton::DressedLepton(const Particle& dlepton)
    : Particle(dlepton)
  {
    setConstituents({{dlepton}}); //< bare lepton is first constituent
  }

  DressedLepton::DressedLepton(const Particle& lepton, const Particles& photons, bool momsum)
    : Particle(lepton.pid(), lepton.momentum())
  {
    setConstituents({{lepton}}); //< bare lepton is first constituent
    addConstituents(photons, momsum);
  }

  void DressedLepton::addPhoton(const Particle& p, bool momsum) {
    if (p.pid() != PID::PHOTON) throw Error("Clustering a non-photon on to a DressedLepton:"+to_string(p.pid()));
    addConstituent(p, momsum);
  }

  const Particle& DressedLepton::bareLepton() const {
    const Particle& l = constituents().front();
    if (!l.isChargedLepton()) throw Error("First constituent of a DressedLepton is not a bare lepton: oops");
    return l;
  }

  //}



  // Separate-FS version
  DressedLeptons::DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                                 double dRmax, const Cut& cut,
                                 bool useDecayPhotons, bool useJetClustering)
    : FinalState(cut),
      _dRmax(dRmax), _fromDecay(useDecayPhotons), _useJetClustering(useJetClustering)
  {
    setName("DressedLeptons");

    // Find photons -- specialising to prompt photons if decay photons are to be vetoed
    IdentifiedFinalState photonfs(photons, PID::PHOTON);
    if (_fromDecay) {
      addProjection(photonfs, "Photons");
    } else {
      addProjection(PromptFinalState(photonfs), "Photons");
    }

    // Find bare leptons
    IdentifiedFinalState leptonfs(bareleptons);
    leptonfs.acceptIdPairs({PID::ELECTRON, PID::MUON, PID::TAU}); //< hmm, no final-state taus, so is this useful?
    addProjection(leptonfs, "Leptons");

    // Set up FJ clustering option
    if (_useJetClustering) {
      MergedFinalState mergedfs(photonfs, leptonfs);
      FastJets leptonjets(mergedfs, FastJets::ANTIKT, dRmax);
      addProjection(leptonjets, "LeptonJets");
    }
  }


  // Single-FS version
  DressedLeptons::DressedLeptons(const FinalState& barefs,
                                 double dRmax, const Cut& cut,
                                 bool useDecayPhotons, bool useJetClustering)
    : DressedLeptons(barefs, barefs, dRmax, cut, useDecayPhotons, useJetClustering)
  {     }




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
            cmp(_fromDecay, other._fromDecay) ||
            cmp(_useJetClustering, other._useJetClustering));
  }


  void DressedLeptons::project(const Event& e) {
    _theParticles.clear();

    // Get bare leptons
    const FinalState& signal = apply<FinalState>(e, "Leptons");
    Particles bareleptons = signal.particles();
    if (bareleptons.empty()) return;

    // Initialise DL collection with bare leptons
    vector<Particle> allClusteredLeptons;
    allClusteredLeptons.reserve(bareleptons.size());

    // If the radius is 0 or negative, don't even attempt to cluster
    if (_useJetClustering) {

      if (_dRmax <= 0) {
        for (const Particle& bl : bareleptons) {
          Particle dl(bl.pid(), bl.momentum());
          dl.setConstituents({bl});
          allClusteredLeptons += dl;
        }
      } else {
        const Jets& lepjets = apply<FastJets>(e, "LeptonJets").jets();
        for (const Jet& lepjet : lepjets) {
          const Particles leps = sortByPt(lepjet.particles(isChargedLepton));
          if (leps.empty()) continue;
          Particles constituents = {leps[0]}; //< note no dressing for subleading leptons
          Particle dl(leps[0].pid(), leps[0].momentum());
          constituents += lepjet.particles(isPhoton);
          dl.setConstituents(constituents);
          allClusteredLeptons += dl;
        }
      }

    } else {

      for (const Particle& bl : bareleptons) {
        Particle dl(bl.pid(), bl.momentum());
        dl.setConstituents({bl});
        allClusteredLeptons += dl;
      }
      if (_dRmax > 0) {
        // Match each photon to its closest charged lepton within the dR cone
        const FinalState& photons = applyProjection<FinalState>(e, "Photons");
        for (const Particle& photon : photons.particles()) {
          // Ignore photon if it's from a hadron/tau decay and we're avoiding those
          /// @todo Can remove via the PromptFinalState conversion above?
          if (!_fromDecay && photon.fromDecay()) continue;
          const FourMomentum& p_P = photon.momentum();
          double dRmin = _dRmax;
          int idx = -1;
          for (size_t i = 0; i < bareleptons.size(); ++i) {
            const Particle& bl = bareleptons[i];
            // Only cluster photons around *charged* signal particles
            if (bl.charge3() == 0) continue;
            // Find the closest lepton
            double dR = deltaR(bl, p_P);
            if (dR < dRmin) {
              dRmin = dR;
              idx = i;
            }
          }
          if (idx > -1) allClusteredLeptons[idx].addConstituent(photon, true);
        }
      }
    }

    // Fill the canonical particles collection with the composite DL Particles
    for (const Particle& lepton : allClusteredLeptons) {
      const bool acc = accept(lepton);
      MSG_TRACE("Clustered lepton " << lepton << " with constituents = " << lepton.constituents() << ", cut-pass = " << boolalpha << acc);
      if (acc) _theParticles.push_back(lepton);
    }
    MSG_DEBUG("#dressed leptons = " << allClusteredLeptons.size() << " -> " << _theParticles.size() << " after cuts");

  }


}
