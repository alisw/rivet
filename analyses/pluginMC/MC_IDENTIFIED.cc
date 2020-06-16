// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  /// @todo Rename as MC_HADRONS
  class MC_IDENTIFIED : public Analysis {
  public:

    /// Constructor
    MC_IDENTIFIED()
      : Analysis("MC_IDENTIFIED")
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState cnfs(Cuts::abseta < 5.0 && Cuts::pT > 500*MeV);
      declare(cnfs, "FS");
      declare(UnstableParticles(Cuts::abseta < 5.0 && Cuts::pT > 500*MeV), "UFS");

      // Histograms
      // @todo Choose E/pT ranged based on input energies... can't do anything about kin. cuts, though

      book(_histStablePIDs  ,"MultsStablePIDs", 3335, -0.5, 3334.5);
      book(_histDecayedPIDs ,"MultsDecayedPIDs", 3335, -0.5, 3334.5);
      book(_histAllPIDs     ,"MultsAllPIDs", 3335, -0.5, 3334.5);

      book(_histEtaPi       ,"EtaPi", 25, 0, 5);
      book(_histEtaK        ,"EtaK", 25, 0, 5);
      book(_histEtaLambda   ,"EtaLambda", 25, 0, 5);
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // Unphysical (debug) plotting of all PIDs in the event, physical or otherwise
      for(ConstGenParticlePtr gp: HepMCUtils::particles(event.genEvent())) {
        _histAllPIDs->fill(abs(gp->pdg_id()), weight);
      }

      // Charged + neutral final state PIDs
      const FinalState& cnfs = apply<FinalState>(event, "FS");
      for (const Particle& p : cnfs.particles()) {
        _histStablePIDs->fill(p.abspid(), weight);
      }

      // Unstable PIDs and identified particle eta spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        _histDecayedPIDs->fill(p.pid(), weight);
        const double eta_abs = p.abseta();
        const PdgId pid = p.abspid(); //if (PID::isMeson(pid) && PID::hasStrange()) {
        if (pid == 211 || pid == 111) _histEtaPi->fill(eta_abs, weight);
        else if (pid == 321 || pid == 130 || pid == 310) _histEtaK->fill(eta_abs, weight);
        else if (pid == 3122) _histEtaLambda->fill(eta_abs, weight);
      }

    }



    /// Finalize
    void finalize() {
      scale(_histStablePIDs, 1/sumOfWeights());
      scale(_histDecayedPIDs, 1/sumOfWeights());
      scale(_histAllPIDs, 1/sumOfWeights());
      scale(_histEtaPi, 1/sumOfWeights());
      scale(_histEtaK, 1/sumOfWeights());
      scale(_histEtaLambda, 1/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histStablePIDs, _histDecayedPIDs, _histAllPIDs;
    Histo1DPtr _histEtaPi, _histEtaK, _histEtaLambda;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_IDENTIFIED);

}
