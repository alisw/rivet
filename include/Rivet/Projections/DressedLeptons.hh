// -*- C++ -*-
#ifndef RIVET_DressedLeptons_HH
#define RIVET_DressedLeptons_HH

#include "Rivet/Tools/Logging.hh"
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Cuts.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// A charged lepton meta-particle created by clustering photons close to the bare lepton
  class DressedLepton : public Particle {
  public:

    DressedLepton(const Particle& lepton) :
      Particle(lepton.pid(), lepton.momentum()),
      _constituentLepton(lepton) {}

    void addPhoton(const Particle& p, bool cluster) {
      _constituentPhotons.push_back(p);
      if (cluster) setMomentum(momentum() + p.momentum());
    }

    const Particle& constituentLepton() const { return _constituentLepton; }
    const Particles& constituentPhotons() const { return _constituentPhotons; }

  private:

    Particles _constituentPhotons;
    Particle _constituentLepton;
  };


  /// @brief Cluster photons from a given FS to all charged particles (typically leptons)
  ///
  /// This stores the original (bare) charged particles and photons as particles()
  /// while the newly created clustered lepton objects are accessible as
  /// dressedLeptons(). The clustering is done by a delta(R) cone around each bare
  /// lepton, with double counting being avoided by only adding a photon to the _closest_
  /// bare lepton if it happens to be within the capture radius of more than one.
  class DressedLeptons : public FinalState {
  public:

    /// @brief Constructor with a general (and optional) Cut argument
    ///
    /// Provide final state projections used to select the photons and bare
    /// leptons (wish we had put the first two args the other way around...),
    /// a clustering delta(R) cone size around each bare lepton, and an optional
    /// cut on the _dressed_ leptons (i.e. the momenta after clustering.)
    /// The final two arguments are rarely used.
    DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                   double dRmax, const Cut& cut=Cuts::open(),
                   bool cluster=true, bool useDecayPhotons=false);

    /// Constructor with a general (and optional) Cut argument
    /// @deprecated Use the version with Cut c before cluster (i.e. with the most common non-default args first)
    DEPRECATED("Use the version with Cut c before cluster")
    DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                   double dRmax, bool cluster=true, const Cut& cut=Cuts::open(),
                   bool useDecayPhotons=false);

    /// Constructor with numerical eta and pT cuts
    /// @deprecated Use the Cut version
    DEPRECATED("Use the Cut version")
    DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                   double dRmax, bool cluster,
                   double etaMin, double etaMax,
                   double pTmin, bool useDecayPhotons=false);


    /// Clone this projection
    virtual const Projection* clone() const {
      return new DressedLeptons(*this);
    }

    /// Retrieve the dressed leptons
    const vector<DressedLepton>& dressedLeptons() const { return _clusteredLeptons; }

    /// Retrieve the dressed leptons (synonym)
    /// @deprecated Use dressedLeptons()
    DEPRECATED("Use dressedLeptons()")
    const vector<DressedLepton>& clusteredLeptons() const { return _clusteredLeptons; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    /// Maximum cone radius to find photons in
    double _dRmax;
    /// Whether to actually add the photon momenta to clusteredLeptons
    bool _cluster;
    /// Whether to include photons from hadron (particularly pi0) decays
    bool _fromDecay;

    /// Container which stores the clustered lepton objects
    vector<DressedLepton> _clusteredLeptons;

  };



}


#endif
