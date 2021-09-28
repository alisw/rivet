// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Isolated photon + jets at 13 TeV
  class ATLAS_2017_I1645627 : public Analysis {
  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1645627);

    // Book histograms and initialise projections before the run
    void init() {
      const FinalState fs;

      // calorimeter particles
      VisibleFinalState visFS(fs);
      VetoedFinalState calo_fs(visFS);
      calo_fs.addVetoPairId(PID::MUON);
      declare(calo_fs, "calo");

      // Voronoi eta-phi tessellation with KT jets, for ambient energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE); // E-scheme used by default;
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::voronoi_area, fastjet::VoronoiAreaSpec(1.0)));
      declare(fj, "KtJetsD05");

      // photon
      PromptFinalState photonfs(Cuts::abspid == PID::PHOTON && Cuts::abseta < 2.37 && Cuts::pT > 125*GeV);
      declare(photonfs, "photons");

      // Jets
      FastJets jetpro(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jetpro, "Jets");

      // Histograms
      book(_h_photon_pt     , 1, 1, 1);
      book(_h_jet_pt        , 2, 1, 1);
      book(_h_phjet_dphi    , 3, 1, 1);
      book(_h_phjet_mass    , 4, 1, 1);
      book(_h_phjet_costheta, 5, 1, 1);

    }


    size_t getEtaBin(double eta) const {
      return binIndex(fabs(eta), _eta_bins_areaoffset);
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the photon
      const Particles& photons = apply<PromptFinalState>(event, "photons").particlesByPt(Cuts::abseta < 1.37 || Cuts::abseta > 1.56);
      if (photons.empty())  vetoEvent;
      const FourMomentum photon = photons[0].momentum();

      // Get the jet
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 100*GeV && Cuts::absrap < 2.37);
      ifilter_discard(jets, deltaRLess(photon, 0.8));
      if (jets.empty())  vetoEvent;
      FourMomentum leadingJet = jets[0].momentum();

      // Compute the jet pT densities
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      FastJets fastjets = apply<FastJets>(event, "KtJetsD05");
      const auto clust_seq_area = fastjets.clusterSeqArea();
      for (const Jet& jet : fastjets.jets()) {
        const double area = clust_seq_area->area(jet); // Implicit call to pseudojet().
        //const double area2 = (clust_seq_area->area_4vector(jet)).perp(); // Area definition used in egammaTruthParticles.
        if (area > 1e-3 && jet.abseta() < _eta_bins_areaoffset.back()) {
          ptDensities.at(getEtaBin(jet.abseta())) += jet.pT()/area;
        }
      }

      // Compute the median event energy density
      vector<double> ptDensity;
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; ++b) {
        ptDensity += ptDensities[b].empty() ? 0 : median(ptDensities[b]);
      }

      // Compute photon isolation with a standard ET cone
      FourMomentum mom_in_EtCone;
      const Particles calo_fs = apply<VetoedFinalState>(event, "calo").particles();
      const double iso_dr = 0.4;
      for (const Particle& p : calo_fs) {
        // Check if it's in the cone of .4
        if (sqrt(2.0*(cosh(p.eta()-photon.eta()) - cos(p.phi()-photon.phi()))) >= iso_dr) continue;
        // Increment sum
        mom_in_EtCone += p.momentum();
      }

      // Remove the photon energy from the isolation
      mom_in_EtCone -= photon;

      // Figure out the correction (area*density)
      const double etcone_area = PI*iso_dr*iso_dr;
      const double correction = ptDensity[getEtaBin(photon.abseta())] * etcone_area;
      // Require photon to be isolated
      if ((mom_in_EtCone.Et()-correction) > (0.0042*photon.pT() + 10*GeV))  vetoEvent;

      // Fill histos
      const double photon_pt = photon.pT()/GeV;
      const double jet_pt = leadingJet.pT()/GeV;
      const double phjet_dphi = deltaPhi(photon, leadingJet);
      const double photon_eta = photon.eta();
      const double jet_y = leadingJet.rapidity();
      _h_photon_pt->fill(photon_pt);
      _h_jet_pt->fill(jet_pt);
      _h_phjet_dphi->fill(phjet_dphi);

      double dy = fabs(jet_y-photon_eta);
      double phjet_costheta = tanh(dy/2.);
      double phjet_mass= (photon+leadingJet).mass()/GeV;
      if (phjet_mass <= 450.)  vetoEvent;
      if (fabs(photon_eta + jet_y) >= 2.37)  vetoEvent;
      if (phjet_costheta >= 0.83)  vetoEvent;
      _h_phjet_costheta->fill(phjet_costheta);
      _h_phjet_mass->fill(phjet_mass);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / sumOfWeights();
      scale(_h_photon_pt, sf);
      scale(_h_jet_pt, sf);
      scale(_h_phjet_dphi, sf);
      scale(_h_phjet_mass, sf);
      scale(_h_phjet_costheta, sf);
    }


  private:

    Histo1DPtr _h_photon_pt, _h_jet_pt;
    Histo1DPtr _h_phjet_dphi, _h_phjet_mass, _h_phjet_costheta;

    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0, 4.0, 5.0};

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1645627);

}
