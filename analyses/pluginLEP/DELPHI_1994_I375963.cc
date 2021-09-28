// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief A study of radiative muon pair events at Z0 energies
  class DELPHI_1994_I375963 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1994_I375963);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // FinalState of prompt photons and muons
      declare(PromptFinalState(Cuts::abspid == PID::MUON && Cuts::energy > 20*GeV), "muons");
      declare(PromptFinalState(Cuts::abspid == PID::PHOTON && Cuts::energy > 2*GeV), "photons");

      // Book histograms
      book(_h["ph_energy"], 1, 1, 1);
      book(_h["ph_angle"], 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      Particles muons = apply<PromptFinalState>(event, "muons").particles();
      Particles photons = apply<PromptFinalState>(event, "photons").particles();
      
      ifilter_select(muons, [](const Particle& muon) {
        double theta = muon.theta()/M_PI * 180.;
        return (theta > 20. && theta < 160.);
      });
      
      ifilter_select(photons, [](const Particle& photon) {
        double theta = photon.theta()/M_PI * 180.;
        double phi = photon.phi()/M_PI * 180.;
        
        bool femc = (theta > 10. && theta < 36.5) || (theta > 143.5 && theta < 170.);
        bool hpc = theta > 43. && theta < 137. && abs(fmod(phi, 15.)) > 1.5 && abs(theta-90.) > 2.;
        
        return (femc || hpc);
      });
            
      if (muons.size() != 2) vetoEvent;
      if (photons.size() == 0) vetoEvent;
      
      for (const Particle& photon : photons) {
        double minAngle = 1000.;
        for (const Particle& muon : muons) {
          double angle = photon.angle(muon) / M_PI * 180.;
          if (angle < minAngle) {
            minAngle = angle;
          }
        }
        if (minAngle > 5.) {
          _h["ph_energy"]->fill(photon.energy()/GeV);
          _h["ph_angle"]->fill(minAngle);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale({_h["ph_energy"], _h["ph_angle"]}, 1./sumW());

    }

    //@}


    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    //@}


  };


  DECLARE_RIVET_PLUGIN(DELPHI_1994_I375963);

}
