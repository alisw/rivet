// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/NeutralFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// Differential cross-section of FSR photons in Z decays
  class CMS_2015_I1346843 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2015_I1346843);

    /// Book histograms and initialise projections before the run
    void init() {

      Cut c_photons = Cuts::pT >= 5*GeV && Cuts::abseta < 2.5 && !(Cuts::absetaIn(1.4, 1.6));
      PromptFinalState photons(Cuts::abspid == PID::PHOTON && c_photons, true, true);
      declare(photons, "PHOTONFS");

      Cut c_muons = Cuts::pT > 9*GeV && Cuts::abseta < 2.4;
      PromptFinalState muons(Cuts::abspid == PID::MUON && c_muons);
      declare(muons, "MUFS");

      book(_h["pho_et"]           ,1, 1, 1);  // photon transverse energy
      book(_h["pho_et_wide"]      ,1, 2, 1);  // photon transverse energy (0.5 < dr < 3.0)
      book(_h["pho_et_close"]     ,1, 3, 1);  // photon transverse energy (0.05 < dr < 0.5)
      book(_h["pho_et_lqt"]       ,1, 4, 1);  // photon transverse energy (q_T < 10)
      book(_h["pho_et_hqt"]       ,1, 5, 1);  // photon transverse energy (q_T > 50)
      book(_h["pho_dr"]           ,2, 1, 1);  // delta_R
      book(_h["pho_dr_lqt"]       ,2, 2, 1);  // delta_R (q_T < 10)
      book(_h["pho_dr_hqt"]       ,2, 3, 1);  // delta_R  (q_T > 50)
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles muons = apply<PromptFinalState>(event, "MUFS").particlesByPt();

      if (muons.size() < 2) vetoEvent;
      if (muons[0].pT() < 31*GeV) vetoEvent;
      if (muons[0].charge()*muons[1].charge() > 0) vetoEvent;
      const double mZ = (muons[0].momentum() + muons[1].momentum()).mass();
      if (!inRange(mZ, 30*GeV, 87*GeV)) vetoEvent;

      const Particles photons = apply<PromptFinalState>(event, "PHOTONFS").particlesByPt();
      // We want the photon with the highest pT that does not come from a decay
      for (const Particle& p : photons) {

        const double dR = std::min(deltaR(p, muons[0]), deltaR(p, muons[1]) );
        if (!inRange(dR, 0.05, 3.0)) continue;

        // Calculate the three-body (mu,mu,gamma) transverse momentum
        const double qT = (muons[0].mom() + muons[1].mom() + p.mom()).pT();

        // Fill the analysis histograms
        _h["pho_et"]->fill(p.pT()/GeV, 1.0);
        _h["pho_dr"]->fill(dR, 1.0);

        _h[(dR <= 0.5 ? "pho_et_close" : "pho_et_wide")]->fill(p.pT()/GeV, 1.0);

        if (qT / GeV < 10.) {
          _h["pho_et_lqt"]->fill(p.pT()/GeV, 1.0);
          _h["pho_dr_lqt"]->fill(dR, 1.0);
        }

        if (qT / GeV > 50.) {
          _h["pho_et_hqt"]->fill(p.pT()/GeV, 1.0);
          _h["pho_dr_hqt"]->fill(dR, 1.0);
        }

        break; // Exit the loop since we found the highest pT lepton already
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection() / sumOfWeights());
    }


  private:

    map<std::string,Histo1DPtr> _h;

  };


  RIVET_DECLARE_PLUGIN(CMS_2015_I1346843);

}
