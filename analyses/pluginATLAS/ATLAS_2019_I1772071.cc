// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Isolated photon + 2 jets at 13 TeV
  class ATLAS_2019_I1772071 : public Analysis {
  public:

    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1772071);

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
      PromptFinalState photonfs(Cuts::abspid == PID::PHOTON && Cuts::abseta < 2.37 && Cuts::pT > 150*GeV);
      declare(photonfs, "photons");

      // Jets
      FastJets jetpro(fs,  FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);       
      declare(jetpro, "Jets");


      vector<string> observables = {"ETGamma", "pTjet", "RapJet","DeltaRapGammaJet", 
                                    "DeltaPhiGammaJet", "DeltaRapJetJet", "DeltaPhiJetJet", 
                                    "MassJetJet", "MassGammaJetJet"};
      vector<string> regions = {"Inclusive","Fragmentation", "Direct"};

      int i=0;
      for (const string& region : regions){ 
        int j = 1;
	      for (const string& name : observables) {
	        book(_h[name+region], 9*i+j, 1, 1);
	        ++j;
	      }
	      ++i;    
      }
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
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 100*GeV && Cuts::absrap < 2.5);
      ifilter_discard(jets, deltaRLess(photon, 0.8));
      if ( jets.size()<2 )  vetoEvent;
      FourMomentum leadingJet = jets[0].momentum();
      FourMomentum subleadingJet = jets[1].momentum();

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
      const double jet_pt1 = leadingJet.pT()/GeV;
      const double jet_pt2 = subleadingJet.pT()/GeV;
      const double jet1_y = leadingJet.rapidity();
      const double jet2_y = subleadingJet.rapidity();
      const double phjet1_dphi = deltaPhi(photon, leadingJet);
      const double phjet2_dphi = deltaPhi(photon, subleadingJet);
      const double phjet1_drap = fabs(photon.eta()-leadingJet.rapidity());
      const double phjet2_drap = fabs(photon.eta()-subleadingJet.rapidity());
      const double jetjet_drap = fabs(leadingJet.rapidity()-subleadingJet.rapidity());
      const FourMomentum jetjet = leadingJet+subleadingJet;
      const double mjetjet = jetjet.mass()/GeV;
      const FourMomentum phjet1 = photon+leadingJet;
      const FourMomentum phjet2 = photon+subleadingJet;
      const FourMomentum phjetjet = photon+leadingJet+subleadingJet;
      const double mphjetjet = phjetjet.mass()/GeV;
      const double jetjet_dphi = deltaPhi(subleadingJet, leadingJet);

      _h["ETGammaInclusive"]->fill(photon_pt);
      _h["pTjetInclusive"]->fill(jet_pt1);
      _h["pTjetInclusive"]->fill(jet_pt2);
      _h["RapJetInclusive"]->fill(fabs(jet1_y));
      _h["RapJetInclusive"]->fill(fabs(jet2_y));
      _h["DeltaRapGammaJetInclusive"]->fill(phjet1_drap);
      _h["DeltaRapGammaJetInclusive"]->fill(phjet2_drap);
      _h["DeltaPhiGammaJetInclusive"]->fill(phjet1_dphi);
      _h["DeltaPhiGammaJetInclusive"]->fill(phjet2_dphi);
      _h["MassJetJetInclusive"]->fill(mjetjet);
      _h["DeltaPhiJetJetInclusive"]->fill(jetjet_dphi);
      _h["DeltaRapJetJetInclusive"]->fill(jetjet_drap);
      _h["MassGammaJetJetInclusive"]->fill(mphjetjet);
      
      if (photon_pt>jet_pt1) {
        _h["ETGammaDirect"]->fill(photon_pt);
        _h["pTjetDirect"]->fill(jet_pt1);
        _h["pTjetDirect"]->fill(jet_pt2);
        _h["RapJetDirect"]->fill(fabs(jet1_y));
        _h["RapJetDirect"]->fill(fabs(jet2_y));
        _h["DeltaRapGammaJetDirect"]->fill(phjet1_drap);
        _h["DeltaRapGammaJetDirect"]->fill(phjet2_drap);
        _h["DeltaPhiGammaJetDirect"]->fill(phjet1_dphi);
        _h["DeltaPhiGammaJetDirect"]->fill(phjet2_dphi);
        _h["MassJetJetDirect"]->fill(mjetjet);
        _h["DeltaPhiJetJetDirect"]->fill(jetjet_dphi);
        _h["DeltaRapJetJetDirect"]->fill(jetjet_drap);
        _h["MassGammaJetJetDirect"]->fill(mphjetjet);

      }
      else if (photon_pt < jet_pt2)	{
        _h["ETGammaFragmentation"]->fill(photon_pt);
        _h["pTjetFragmentation"]->fill(jet_pt1);
        _h["pTjetFragmentation"]->fill(jet_pt2);
        _h["RapJetFragmentation"]->fill(fabs(jet1_y));
        _h["RapJetFragmentation"]->fill(fabs(jet2_y));
        _h["DeltaRapGammaJetFragmentation"]->fill(phjet1_drap);
        _h["DeltaRapGammaJetFragmentation"]->fill(phjet2_drap);
        _h["DeltaPhiGammaJetFragmentation"]->fill(phjet1_dphi);
        _h["DeltaPhiGammaJetFragmentation"]->fill(phjet2_dphi);
        _h["MassJetJetFragmentation"]->fill(mjetjet);
        _h["DeltaPhiJetJetFragmentation"]->fill(jetjet_dphi);
        _h["DeltaRapJetJetFragmentation"]->fill(jetjet_drap);
        _h["MassGammaJetJetFragmentation"]->fill(mphjetjet);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / picobarn / sumOfWeights();
	    scale(_h, sf);
    }


  private:

    map<string,Histo1DPtr> _h;
    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0, 4.0, 5.0};

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1772071);

}
