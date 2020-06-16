// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  


  /// @brief MC validation analysis for single photon events
  class MC_PHOTONINC : public Analysis {
  public:

    /// Default constructor
    MC_PHOTONINC()
      : Analysis("MC_PHOTONINC")
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

      // FS for isolation excludes the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "JetFS");

      book(_h_photon_pT ,"photon_pT", logspace(50, 30.0, 0.5*(sqrtS()>0.?sqrtS():14000.)));
      book(_h_photon_pT_lin ,"photon_pT_lin", 50, 0.0, 70.0);
      book(_h_photon_y ,"photon_y", 50, -5.0, 5.0);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      // Get the photon
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

      _h_photon_pT->fill(photon.pT(),weight);
      _h_photon_pT_lin->fill(photon.pT(),weight);
      _h_photon_y->fill(photon.rapidity(),weight);
    }


    // Finalize
    void finalize() {
      scale(_h_photon_pT, crossSectionPerEvent());
      scale(_h_photon_pT_lin, crossSectionPerEvent());
      scale(_h_photon_y, crossSectionPerEvent());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_photon_pT;
    Histo1DPtr _h_photon_pT_lin;
    Histo1DPtr _h_photon_y;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_PHOTONINC);

}
