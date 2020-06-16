// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  


  /// @brief MC validation analysis for photon + jets events
  class MC_PHOTONJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_PHOTONJETS()
      : MC_JetAnalysis("MC_PHOTONJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      // General FS
      FinalState fs((Cuts::etaIn(-5.0, 5.0)));
      declare(fs, "FS");

      // Get leading photon
      LeadingParticlesFinalState photonfs(FinalState((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  30.0*GeV)));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // FS for jets excludes the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "JetFS");
      FastJets jetpro(vfs, FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      book(_h_photon_jet1_deta ,"photon_jet1_deta", 50, -5.0, 5.0);
      book(_h_photon_jet1_dphi ,"photon_jet1_dphi", 20, 0.0, M_PI);
      book(_h_photon_jet1_dR ,"photon_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }


    /// Do the analysis
    void analyze(const Event& e) {
      // Get the photon
      /// @todo share IsolatedPhoton projection between all MC_*PHOTON* analyses
      const Particles photons = apply<FinalState>(e, "LeadingPhoton").particles();
      if (photons.size() != 1) {
        vetoEvent;
      }
      const FourMomentum photon = photons.front().momentum();

      // Get all charged particles
      const FinalState& fs = apply<FinalState>(e, "JetFS");
      if (fs.empty()) {
        vetoEvent;
      }

      // Passed cuts, so get the weight
      const double weight = 1.0;

      // Isolate photon by ensuring that a 0.4 cone around it contains less than 7% of the photon's energy
      const double egamma = photon.E();
      double econe = 0.0;
      for (const Particle& p : fs.particles()) {
        if (deltaR(photon, p.momentum()) < 0.4) {
          econe += p.E();
          // Veto as soon as E_cone gets larger
          if (econe/egamma > 0.07) {
            vetoEvent;
          }
        }
      }

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(_jetptcut);
      if (jets.size()>0) {
        _h_photon_jet1_deta->fill(photon.eta()-jets[0].eta(), weight);
        _h_photon_jet1_dphi->fill(mapAngle0ToPi(photon.phi()-jets[0].phi()), weight);
        _h_photon_jet1_dR->fill(deltaR(photon, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    // Finalize
    void finalize() {
      scale(_h_photon_jet1_deta, crossSectionPerEvent());
      scale(_h_photon_jet1_dphi, crossSectionPerEvent());
      scale(_h_photon_jet1_dR, crossSectionPerEvent());

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_photon_jet1_deta;
    Histo1DPtr _h_photon_jet1_dphi;
    Histo1DPtr _h_photon_jet1_dR;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_PHOTONJETS);

}
