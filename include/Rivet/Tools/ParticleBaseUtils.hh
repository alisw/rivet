#ifndef RIVET_PARTICLEBASEUTILS_HH
#define RIVET_PARTICLEBASEUTILS_HH

#include "Rivet/ParticleBase.hh"

namespace Rivet {


  /// @name ParticleBase classifying functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to any() or all() e.g. any(jets, DeltaRLess(electron, 0.4))
  //@{

  /// Base type for Particle -> bool functors
  struct BoolParticleBaseFunctor {
    virtual bool operator()(const ParticleBase& p) const = 0;
  };

  /// Transverse momentum greater-than functor
  struct pTGtr : public BoolParticleBaseFunctor {
    pTGtr(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() > ptcut; }
    double ptcut;
  };
  using PtGtr = pTGtr;

  /// Transverse momentum less-than functor
  struct pTLess : public BoolParticleBaseFunctor {
    pTLess(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() < ptcut; }
    double ptcut;
  };
  using PtLess = pTLess;


  /// Pseudorapidity greater-than functor
  struct etaGtr : public BoolParticleBaseFunctor {
    etaGtr(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() > etacut; }
    double etacut;
  };
  using EtaGtr = etaGtr;

  /// Pseudorapidity momentum less-than functor
  struct etaLess : public BoolParticleBaseFunctor {
    etaLess(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() < etacut; }
    double etacut;
  };
  using EtaLess = etaLess;

  /// Abs pseudorapidity greater-than functor
  struct absetaGtr : public BoolParticleBaseFunctor {
    absetaGtr(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() > absetacut; }
    double absetacut;
  };
  using AbsEtaGtr = etaGtr;

  /// Abs pseudorapidity momentum less-than functor
  struct absetaLess : public BoolParticleBaseFunctor {
    absetaLess(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() < absetacut; }
    double absetacut;
  };
  using AbsEtaLess = absetaLess;


  /// Rapidity greater-than functor
  struct rapGtr : public BoolParticleBaseFunctor {
    rapGtr(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() > rapcut; }
    double rapcut;
  };
  using RapGtr = rapGtr;

  /// Rapidity momentum less-than functor
  struct rapLess : public BoolParticleBaseFunctor {
    rapLess(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() < rapcut; }
    double rapcut;
  };
  using RapLess = rapLess;

  /// Abs rapidity greater-than functor
  struct absrapGtr : public BoolParticleBaseFunctor {
    absrapGtr(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() > absrapcut; }
    double absrapcut;
  };
  using AbsRapGtr = absrapGtr;

  /// Abs rapidity momentum less-than functor
  struct absrapLess : public BoolParticleBaseFunctor {
    absrapLess(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() < absrapcut; }
    double absrapcut;
  };
  using AbsRapLess = absrapLess;


  /// Delta R (with respect to another 4-momentum, @a vec) greater-than functor
  struct deltaRGtr : public BoolParticleBaseFunctor {
    deltaRGtr(const FourMomentum& vec, double dr) : refvec(vec), drcut(dr) { }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec) > drcut; }
    FourMomentum refvec;
    double drcut;
  };
  using DeltaRGtr = deltaRGtr;

  /// Delta R (with respect to another 4-momentum, @a vec) less-than functor
  struct deltaRLess : public BoolParticleBaseFunctor {
    deltaRLess(const FourMomentum& vec, double dr) : refvec(vec), drcut(dr) { }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec) < drcut; }
    FourMomentum refvec;
    double drcut;
  };
  using DeltaRLess = deltaRLess;

  //@}


  /// @name Non-PID particle properties, via unbound functions
  /// @todo Move to FourMomentum functions
  //@{

  /// Unbound function access to momentum
  inline FourMomentum mom(const ParticleBase& p) { return p.mom(); }

  /// Unbound function access to p3
  inline Vector3 p3(const ParticleBase& p) { return p.p3(); }

  /// Unbound function access to p
  inline double p(const ParticleBase& p) { return p.p(); }

  /// Unbound function access to pT
  inline double pT(const ParticleBase& p) { return p.pT(); }

  /// Unbound function access to eta
  inline double eta(const ParticleBase& p) { return p.eta(); }

  /// Unbound function access to abseta
  inline double abseta(const ParticleBase& p) { return p.abseta(); }

  /// Unbound function access to rapidity
  inline double rap(const ParticleBase& p) { return p.rap(); }

  /// Unbound function access to abs rapidity
  inline double absrap(const ParticleBase& p) { return p.absrap(); }

  //@}


}

#endif
