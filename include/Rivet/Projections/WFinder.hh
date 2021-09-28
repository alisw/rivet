// -*- C++ -*-
#ifndef RIVET_WFinder_HH
#define RIVET_WFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying W
  ///
  /// Chain together different projections as convenience for finding one W
  /// from one lepton and the missing E 4-vector in the final state, including photon clustering.
  class WFinder : public ParticleFinder {
  public:

    enum class ChargedLeptons { PROMPT, ALL };
    enum class ClusterPhotons { NONE, NODECAY, ALL };
    enum class AddPhotons { NO, YES };
    enum class MassWindow { M, MT };


    /// @name Constructors
    //@{

    /// @brief Constructor taking cuts object
    ///
    /// @param inputfs  Input final state
    /// @param leptoncuts  Charged lepton cuts
    /// @param pid  Type of the charged lepton
    /// @param minmass,maxmass  (Transverse) mass window
    /// @param missingET  Minimal amount of missing ET (neutrinos) required
    /// @param dRmax  Maximum dR of photons around charged lepton to take into account
    ///  for W reconstruction (only relevant if one of the following are true)
    /// @param chLeptons  Only use prompt charged leptons, or any charged leptons?
    /// @param clusterPhotons  Whether such photons are supposed to be
    ///  clustered to the lepton object and thus W mom
    /// @param trackPhotons  Whether such photons should be added to _theParticles
    /// @param masstype  Whether mass window should be applied using m or mT
    /// @param masstarget  The expected (transverse) mass value, if resolving ambiguities
    ///
    /// @todo Revisit AddPhotons::NO as default?
    WFinder(const FinalState& inputfs,
            const Cut& leptoncuts,
            PdgId pid,
            double minmass, double maxmass,
            double missingET,
            double dRmax=0.1,
            ChargedLeptons chLeptons=ChargedLeptons::PROMPT,
            ClusterPhotons clusterPhotons=ClusterPhotons::NODECAY,
            AddPhotons trackPhotons=AddPhotons::NO,
            MassWindow masstype=MassWindow::M,
            double masstarget=80.4*GeV);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(WFinder);

    //@}


    /// @brief Access to the found bosons, equivalent to constituents()
    /// @note Currently either 0 or 1 boson can be found.
    const Particles& bosons() const { return particles(); }
    /// Access to the found boson (assuming it exists)
    /// @todo C++17 std::optional...
    const Particle& boson() const { return particles().front(); }


    /// @brief Access to the Ws' constituent clustered leptons
    /// @note Either size 0 if no boson was found or 1 if one boson was found
    const Particles& constituentLeptons() const { return _leptons; }
    /// brief Access to the W's constituent clustered lepton (assuming it exists)
    /// @todo C++17 std::optional...
    const Particle& constituentLepton() const { return _leptons.front(); }


    /// Access to the Ws' constituent neutrinos
    ///
    /// @note Either size 0 if no boson was found or 1 if one boson was found
    /// @note The neutrino can't be perfecly reconstructed -- this is a pseudo-nu from the MET.
    const Particles& constituentNeutrinos() const { return _neutrinos; }
    /// Access to the W's constituent neutrino (assuming it exists)
    /// @note The neutrino can't be perfecly reconstructed -- this is a pseudo-nu from the MET.
    const Particle& constituentNeutrino() const { return _neutrinos.front(); }


    /// Access to the particles other than the W leptons and clustered photons
    ///
    /// Useful for e.g. input to a jet finder
    const VetoedFinalState& remainingFinalState() const;

    /// Access to the missing momentum projection used to find the "neutrino"
    const MissingMomentum& missingMom() const;

    /// @brief Calculate the transverse mass of the W, from the charged lepton and neutrino
    ///
    /// Defined as sqrt(2 pT_l pT_nu (1.0 - cos(dphi_lnu))). Return -1 if no boson found.
    double mT() const {
      if (bosons().empty()) return -1;
      return Rivet::mT(constituentLepton().mom(), constituentNeutrino().mom());
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    CmpState compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() { _theParticles.clear(); }


  private:

    /// (Transverse) mass cuts
    double _minmass, _maxmass, _masstarget;

    /// Use transverse or complete mass?
    bool _useTransverseMass;

    /// Missing ET cut
    double _etMissMin;

    /// Switch for tracking of photons (whether to include them in the W particle)
    /// This is relevant when the clustered photons need to be excluded from e.g. a jet finder
    AddPhotons _trackPhotons;

    /// Charged lepton flavour
    PdgId _pid;

    /// Result caches. Will be filled by project()
    Particles _leptons, _neutrinos;

  };


}


#endif
