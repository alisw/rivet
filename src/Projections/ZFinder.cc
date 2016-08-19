// -*- C++ -*-
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  ZFinder::ZFinder(const FinalState& inputfs,
		   const Cut & fsCut,
                   PdgId pid,
                   double minmass, double maxmass,
                   double dRmax,
                   ClusterPhotons clusterPhotons,
                   PhotonTracking trackPhotons,
                   double masstarget)
  {
    setName("ZFinder");

    _minmass = minmass;
    _maxmass = maxmass;
    _masstarget = masstarget;
    _pid = pid;
    _trackPhotons = trackPhotons;

    IdentifiedFinalState bareleptons(inputfs);
    bareleptons.acceptIdPair(pid);
    const bool doClustering = (clusterPhotons != NOCLUSTER);
    const bool useDecayPhotons = (clusterPhotons == CLUSTERALL);
    DressedLeptons leptons(inputfs, bareleptons, dRmax, fsCut, doClustering, useDecayPhotons);
    addProjection(leptons, "DressedLeptons");

    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    addProjection(remainingFS, "RFS");
  }


  /////////////////////////////////////////////////////


  const VetoedFinalState& ZFinder::remainingFinalState() const {
    return getProjection<VetoedFinalState>("RFS");
  }


  int ZFinder::compare(const Projection& p) const {
    PCmp LCcmp = mkNamedPCmp(p, "DressedLeptons");
    if (LCcmp != EQUIVALENT) return LCcmp;

    const ZFinder& other = dynamic_cast<const ZFinder&>(p);
    return (cmp(_minmass, other._minmass) || cmp(_maxmass, other._maxmass) ||
            cmp(_pid, other._pid) || cmp(_trackPhotons, other._trackPhotons));
  }


  void ZFinder::project(const Event& e) {
    clear();

    const DressedLeptons& leptons = applyProjection<DressedLeptons>(e, "DressedLeptons");

    InvMassFinalState imfs(std::make_pair(_pid, -_pid), _minmass, _maxmass, _masstarget);
    Particles tmp;
    tmp.insert(tmp.end(), leptons.dressedLeptons().begin(), leptons.dressedLeptons().end());
    imfs.calc(tmp);

    if (imfs.particlePairs().size() < 1) return;
    ParticlePair Zconstituents(imfs.particlePairs()[0]);
    Particle l1(Zconstituents.first), l2(Zconstituents.second);
    if (threeCharge(l1) > 0) {
      _constituents.push_back(l1);
      _constituents.push_back(l2);
    } else {
      _constituents.push_back(l2);
      _constituents.push_back(l1);
    }
    FourMomentum pZ = l1.momentum() + l2.momentum();
    assert(threeCharge(l1) + threeCharge(l2) == 0);

    stringstream msg;
    msg << "Z " << pZ << " reconstructed from: \n"
        << "   " << l1.momentum() << " " << l1.pid() << "\n"
        << " + " << l2.momentum() << " " << l2.pid();
    MSG_DEBUG(msg.str());
    _bosons.push_back(Particle(PID::ZBOSON, pZ));

    // Keep the DressedLeptons found by the ZFinder
    foreach (DressedLepton l, leptons.dressedLeptons()) {
      _allLeptons.push_back(l);
    }

    // Find the DressedLeptons which survived the IMFS cut such that we can
    // extract their original particles
    foreach (const Particle& p, _constituents) {
      foreach (const DressedLepton& l, leptons.dressedLeptons()) {
        if (p.pid() == l.pid() && p.momentum() == l.momentum()) {
          _theParticles.push_back(l.constituentLepton());
          if (_trackPhotons) {
            _theParticles.insert(_theParticles.end(), l.constituentPhotons().begin(), l.constituentPhotons().end());
          }
        }
      }
    }
  }


}
