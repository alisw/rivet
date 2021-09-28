// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief D0 inclusive isolated photon cross-section vs. \f$ p_\perp(gamma) \f$.
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  class D0_2006_S6438750 : public Analysis {

  public:

    /// @name Constructors etc.
    //@{

    /// Default constructor.
    D0_2006_S6438750()
      : Analysis("D0_2006_S6438750")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    void init() {
      // General FS for photon isolation
      FinalState fs;
      declare(fs, "AllFS");

      // Get leading photon
      LeadingParticlesFinalState photonfs(FinalState((Cuts::etaIn(-0.9, 0.9) && Cuts::pT >=  23.0*GeV)));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // Book histograms
      book(_h_pTgamma ,1, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // Get the photon
      const FinalState& photonfs = apply<FinalState>(event, "LeadingPhoton");
      if (photonfs.particles().size() != 1) {
        vetoEvent;
      }
      const FourMomentum photon = photonfs.particles().front().momentum();

      // Isolate photon by ensuring that a 0.4 cone around it contains less than 10% of the photon's energy
      double E_P   = photon.E();
      double eta_P = photon.eta();
      double phi_P = photon.phi();
      double econe = 0.0;
      for (const Particle& p : apply<FinalState>(event, "AllFS").particles()) {
        if (deltaR(eta_P, phi_P,
                   p.eta(), p.phi()) < 0.4) {
          econe += p.E();
          if (econe/E_P > 1.1) {
            vetoEvent;
          }
        }
      }

      // Fill histo
      _h_pTgamma->fill(photon.pT());
    }



    // Finalize
    void finalize() {
      const double lumi_gen = sumOfWeights()/crossSection();
      // Divide by effective lumi, plus rapidity bin width of 1.8
      scale(_h_pTgamma, 1/lumi_gen * 1/1.8);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_pTgamma;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2006_S6438750);

}
