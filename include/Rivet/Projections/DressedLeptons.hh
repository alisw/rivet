// -*- C++ -*-
#ifndef RIVET_DressedLeptons_HH
#define RIVET_DressedLeptons_HH

#include "Rivet/Projection.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Config/RivetCommon.hh"

namespace Rivet {


  /// @brief A charged lepton meta-particle created by clustering photons close to the bare lepton
  ///
  /// @todo Remove completely -- it's unnecessary and too confusing (esp. between copying & aggregating)
  /// @deprecated Just use Particle.constituents() now.
  class DressedLepton : public Particle {
  public:

    /// Copy constructor (from Particle)
    DressedLepton(const Particle& dlepton);

    /// Components constructor
    /// @note This is not a copy constructor, hence the explicit second argument even if empty
    DressedLepton(const Particle& lepton, const Particles& photons, bool momsum=true);

    /// Add a photon to the dressed lepton
    /// @todo Deprecate and override add/setConstituents instead?
    void addPhoton(const Particle& p, bool momsum=true);

    /// Retrieve the bare lepton
    const Particle& bareLepton() const;
    /// Retrieve the bare lepton (alias)
    /// @deprecated Prefer the more physicsy bareLepton()
    const Particle& constituentLepton() const { return bareLepton(); }

    /// Retrieve the clustered photons
    const Particles photons() const { return slice(constituents(), 1); }
    /// Retrieve the clustered photons (alias)
    /// @deprecated Prefer the shorter photons()
    const Particles constituentPhotons() const { return photons(); }

  };


  /// @brief Cluster photons from a given FS to all charged particles (typically leptons)
  ///
  /// The clustering is done by a delta(R) cone around each bare lepton or by
  /// jet clustering. In both modes, double counting is avoided: for the dR
  /// clustering, a photon is only added to the _closest_ bare lepton if it
  /// happens to be within the capture radius of more than one; for the jet
  /// clustering, only the bare lepton with the highest pT is retained if more
  /// than one is clustered into a jet.
  ///
  /// @note The particles() and dressedLeptons() methods both return the
  /// composite clustered-lepton objects, just with a few extra helper methods
  /// on the special DressedLepton type returned by the latter. The constituent
  /// bare leptons and photons are returned by rawParticles() (inherited from
  /// ParticleFinder)
  ///
  class DressedLeptons : public FinalState {
  public:

    /// @brief Constructor with a single input FinalState (used for both photons and bare leptons)
    ///
    /// Provide a single final state projection used to select the photons and
    /// bare leptons, a photon-clustering delta(R) cone size around each bare
    /// lepton, and an optional cut on the _dressed_ leptons (i.e. the momenta
    /// and PID after clustering).  The final arguments control whether
    /// non-prompt photons are to be included, and whether the matching of
    /// photons to leptons is to be done via dR matching to the bare lepton or
    /// by a jet clustering algorithm.  Set the clustering radius to 0 or
    /// negative to disable clustering.
    DressedLeptons(const FinalState& allfs,
                   double dRmax, const Cut& cut=Cuts::open(),
                   bool useDecayPhotons=false,
                   bool useJetClustering=false);

    /// @brief Constructor with default input FinalState
    ///
    /// DressedLepton construction from a default-constructed FinalState.
    /// Provide a photon-clustering delta(R) cone size around each bare lepton,
    /// and an optional cut on the _dressed_ leptons (i.e. the momenta and PID
    /// after clustering).  The final arguments control whether non-prompt
    /// photons are to be included, and whether the matching of photons to
    /// leptons is to be done via dR matching to the bare lepton or by a jet
    /// clustering algorithm.  Set the clustering radius to 0 or negative to
    /// disable clustering.
    DressedLeptons(double dRmax, const Cut& cut=Cuts::open(),
                   bool useDecayPhotons=false,
                   bool useJetClustering=false)
      : DressedLeptons(FinalState(), dRmax, cut, useDecayPhotons, useJetClustering)
    {   }

    /// @brief Constructor with distinct photon and lepton finders
    ///
    /// Provide final state projections used to select the photons and bare
    /// leptons, a clustering delta(R) cone size around each bare lepton, and an
    /// optional cut on the _dressed_ leptons (i.e. the momenta and PID after
    /// clustering.)  The final arguments control whether non-prompt photons are
    /// to be included, and whether the matching of photons to leptons is to be
    /// done via dR matching to the bare lepton or by a jet clustering
    /// algorithm.  Set the clustering radius to 0 or negative to disable
    /// clustering.
    ///
    /// @note Wish we had put the first two args the other way around...
    ///
    /// @todo Convert second arg to a general ParticleFinder rather than an FS, to
    /// allow clustering on to unstables, e.g. taus via TauFinder.
    DressedLeptons(const FinalState& photons, const FinalState& bareleptons,
                   double dRmax, const Cut& cut=Cuts::open(),
                   bool useDecayPhotons=false,
                   bool useJetClustering=false);


    /// Clone this projection
    DEFAULT_RIVET_PROJ_CLONE(DressedLeptons);


    /// @brief Retrieve the dressed leptons
    ///
    /// @note Like particles() but with helper functions
    vector<DressedLepton> dressedLeptons() const {
      vector<DressedLepton> rtn;
      for (const Particle& p : particles(cmpMomByPt))
        rtn += DressedLepton(p);  //static_cast<const DressedLepton>(p);
      return rtn;
    }

    /// @brief Retrieve the dressed leptons ordered by supplied sorting functor
    ///
    /// @note Like particles() but with helper functions
    vector<DressedLepton> dressedLeptons(const ParticleSorter& sorter) const {
      vector<DressedLepton> rtn;
      for (const Particle& p : particles(sorter))
        rtn += DressedLepton(p);  //static_cast<const DressedLepton>(p);
      return rtn;
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  private:

    /// Maximum cone radius to find photons in
    double _dRmax;

    /// Whether to include photons from hadron (particularly pi0) and hadronic tau decays
    bool _fromDecay;

    /// Whether to use a jet clustering algorithm rather than nearest-lepton association
    bool _useJetClustering;


  };



}


#endif
