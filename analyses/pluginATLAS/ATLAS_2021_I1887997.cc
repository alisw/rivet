// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Isolated diphoton + X differential cross-sections with full run-2
  class ATLAS_2021_I1887997 : public Analysis {
  public:

    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2021_I1887997);


    // Book histograms and initialise projections before the run
    void init() {

      // Calorimeter particles for photon isolation
      VisibleFinalState visFS;
      VetoedFinalState calo_fs(visFS);
      calo_fs.addVetoPairId(PID::MUON);
      declare(calo_fs, "calo_fs");

      // Photons
      declare(PromptFinalState(Cuts::abspid == PID::PHOTON), "Photons");

      // Jets for UE subtraction with jet-area method
      FastJets fj(FinalState(), FastJets::KT, 0.5, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec(0.9)));
      declare(fj, "KtJetsD05");

      // Histograms
      book(_xs, "yy_xs");
      _observables = {"ph1_pt","ph2_pt","yy_cosTS","yy_m",
                      "yy_phiStar","yy_piMDphi","yy_pT","yy_pTt"};
      for (auto name : _observables) {
        book(_h[name], name);
      }
    }

    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Require at least 2 prompt photons in final state
      Particles photons = apply<PromptFinalState>(event, "Photons").particlesByPt();
      if (photons.size() < 2) vetoEvent;
      photons.resize(2);
      const FourMomentum ph1 = photons[0];
      const FourMomentum ph2 = photons[1];

      // Leading photon should have pT > 40 GeV, subleading > 30 GeV
      const double ph1_pt = ph1.pT();
      const double ph2_pt = ph2.pT();
      if (ph1_pt < 40.*GeV || ph2_pt < 30.*GeV) vetoEvent;

      // Apply photon eta cuts
      ifilter_select(photons, (Cuts::abseta < 2.37) && ( (Cuts::abseta <= 1.37) || (Cuts::abseta >= 1.52) ));
      if (photons.size() < 2) vetoEvent;

      // Require the two photons to be separated in dR
      if (deltaR(ph1,ph2) < 0.4) vetoEvent;

      // Get UE pt densities rho for subtraction later
      const vector<double> eta_bins = {0.0, 1.5, 3.0};
      vector<double> rho(eta_bins.size()-1, 0.0);
      FastJets ktjets = applyProjection<FastJets>(event, "KtJetsD05");
      for (size_t ieta = 0; ieta < eta_bins.size()-1; ++ieta) {
        fastjet::Selector fjselector(fastjet::SelectorAbsRapRange(eta_bins[ieta], eta_bins[ieta+1]));
        double sigma, area;
        ktjets.clusterSeqArea()->get_median_rho_and_sigma(fjselector, true,
                                                          rho[ieta], sigma, area);
      }

      // Loop over photons and require isolation
      const double isoRCone=0.2;
      for (const Particle& photon : photons) {
        // Compute calo isolation via particles within a cone around the photon
        const Particles fs = apply<VetoedFinalState>(event, "calo_fs").particles();
        FourMomentum mom_in_EtCone;
        for (const Particle& p : fs) {
          // Reject if not in cone
          if (deltaR(photon.momentum(), p.momentum()) > isoRCone)  continue;
          // Sum momentum
          mom_in_EtCone += p.momentum();
        }
        // subtract core photon
        mom_in_EtCone -= photon.momentum();
        // UE subtraction energy
        double UEpT = M_PI*sqr(isoRCone) * rho[binIndex(fabs(photon.eta()), eta_bins)];
        // Use photon if energy in isolation cone is low enough
        if (mom_in_EtCone.Et() - UEpT > 0.09*photon.momentum().pT()) vetoEvent;
      }

      map<string, double> obs;
      obs["ph1_pt"] = ph1_pt;
      obs["ph2_pt"] = ph2_pt;
      const FourMomentum yy = ph1 + ph2;
      obs["yy_m"] = yy.mass();
      obs["yy_pT"] = yy.pT();
      obs["yy_piMDphi"] = PI-mapAngle0ToPi(ph1.phi() - ph2.phi());

      obs["yy_cosTS"] = fabs(sinh(( ph1.eta() - ph2.eta() ))*2.0*ph1_pt*ph2_pt/sqrt(sqr(obs["yy_m"])+sqr(obs["yy_pT"]))/obs["yy_m"]); // Collins Soper frame
      const double yy_cosTSLab = fabs(tanh(( ph1.eta() - ph2.eta() ) / 2.)); // Lab frame
      const double sinthetastar_ = sqrt(1. - pow(yy_cosTSLab, 2));
      obs["yy_phiStar"] = tan(0.5 * obs["yy_piMDphi"]) * sinthetastar_;

      // a_t
      const Vector3 t_hat(ph1.x()-ph2.x(), ph1.y()-ph2.y(), 0.);
      const double factor = t_hat.mod();
      const Vector3 t_hatx(t_hat.x()/factor, t_hat.y()/factor, t_hat.z()/factor);
      const Vector3 At(ph1.x()+ph2.x(), ph1.y()+ph2.y(), 0.);
      // Compute a_t transverse component with respect to t_hat
      obs["yy_pTt"] = At.cross(t_hatx).mod();

      // Fill fiducial cross section
      _xs->fill();

      // Fill histograms
      for (auto name : _observables) {
        _h[name]->fill(obs[name]);
      }
    }


    // Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / (picobarn * sumOfWeights());
      scale(_xs, sf);
      scale(_h, sf);
    }


  private:

    CounterPtr _xs;
    map<string, Histo1DPtr> _h;
    vector<string> _observables;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2021_I1887997);

}
