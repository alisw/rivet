// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  /// @brief CDF inclusive isolated prompt photon cross-section
  class CDF_2009_S8436959 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2009_S8436959);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      declare(fs, "FS");

      LeadingParticlesFinalState photonfs(FinalState((Cuts::etaIn(-1.0, 1.0) && Cuts::pT >=  30.0*GeV)));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      book(_h_Et_photon ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles fs = apply<FinalState>(event, "FS").particles();
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size()!=1) {
        vetoEvent;
      }
      FourMomentum leadingPhoton = photons[0].momentum();
      double eta_P = leadingPhoton.eta();
      double phi_P = leadingPhoton.phi();
      FourMomentum mom_in_cone;
      for (const Particle& p : fs) {
        if (deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.4) {
            mom_in_cone += p.momentum();
        }
      }
      if ( (mom_in_cone.Et() - leadingPhoton.Et()) > 2.0*GeV) {
        vetoEvent;
      }
      _h_Et_photon->fill(leadingPhoton.Et());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Et_photon, crossSection()/sumOfWeights()/2.0);
    }

    //@}


  private:

    /// Histogram
    Histo1DPtr _h_Et_photon;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2009_S8436959, CDF_2009_I834437);

}
