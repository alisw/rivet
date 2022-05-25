// -*- C++ -*-
#ifndef RIVET_ZFinder_HH
#define RIVET_ZFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// Chain together different projections as convenience for finding Z's
  /// from two leptons in the final state, including photon clustering.
  ///
  /// @todo Alias then rename as Dileptons
  class ZFinder : public ParticleFinder {
  public:

    enum class ChargedLeptons { PROMPT, ALL };
    enum class ClusterPhotons { NONE, NODECAY, ALL };
    enum class AddPhotons { NO, YES };

    /// @name Constructors
    //@{

    /// @brief Constructor taking cuts object
    ///
    /// @param inputfs Input final state
    /// @param cuts  Lepton cuts
    /// @param pid  Type of the leptons
    /// @param minmass,maxmass  Dilepton mass window
    /// @param dRmax  Maximum dR of photons around leptons to take into account
    ///  for Z reconstruction (only relevant if one of the following are true)
    /// @param chLeptons  The type of charged leptons considered
    /// @param clusterPhotons  Whether such photons are supposed to be
    ///  clustered to the lepton objects and thus Z mom
    /// @param trackPhotons  Whether such photons should be considered constituent particles
    /// @param masstarget  The expected (transverse) mass value, if resolving ambiguities
    ZFinder(const FinalState& inputfs,
            const Cut& cuts,
            PdgId pid,
            double minmass, double maxmass,
            double dRmax=0.1,
            ChargedLeptons chLeptons=ChargedLeptons::PROMPT,
            ClusterPhotons clusterPhotons=ClusterPhotons::NODECAY,
            AddPhotons trackPhotons=AddPhotons::NO,
            double masstarget=91.2*GeV);

    /// Backward-compatible constructor with implicit chLeptons mode = PROMPTCHLEPTONS
    /// @deprecated Remove this and always use the constructor with chLeptons argument.
    ZFinder(const FinalState& inputfs,
            const Cut& cuts,
            PdgId pid,
            double minmass, double maxmass,
            double dRmax,
            ClusterPhotons clusterPhotons,
            AddPhotons trackPhotons=AddPhotons::NO,
            double masstarget=91.2*GeV)
      : ZFinder(inputfs, cuts, pid, minmass, maxmass,
                dRmax, ChargedLeptons::PROMPT, clusterPhotons, trackPhotons, masstarget)
    {   }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(ZFinder);

    //@}


    /// Access to the found bosons
    ///
    /// @note Currently either 0 or 1 boson can be found.
    const Particles& bosons() const { return particles(); }
    /// Access to the found boson (assuming it exists).
    const Particle& boson() const { return bosons().front(); }


    /// Access to the Z constituent clustered leptons
    ///
    /// For example, to make more fine-grained cuts on the clustered leptons.
    /// The positive charge constituent is first in the list (if not empty), and
    /// the negative one second.
    const Particles & constituentLeptons() const;
    const Particles & constituents() const { return constituentLeptons(); }

    /// Access to the particles other than the Z leptons and clustered photons
    ///
    /// Useful for e.g. input to a jet finder
    const VetoedFinalState& remainingFinalState() const;


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() { _theParticles.clear(); }


  private:

    /// Mass cuts to apply to clustered leptons (cf. InvMassFinalState)
    double _minmass, _maxmass, _masstarget;

    /// Switch for tracking of photons (whether to include them in the Z particle)
    /// This is relevant when the clustered photons need to be excluded from e.g. a jet finder
    AddPhotons _trackPhotons;

    /// Lepton flavour
    PdgId _pid;

  };


}

#endif
