// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Inclusive isolated prompt photon analysis with 2011 LHC data
  class ATLAS_2013_I1263495 : public Analysis {
  public:

    /// Constructor
    ATLAS_2013_I1263495()
      : Analysis("ATLAS_2013_I1263495")
    {
      _eta_bins.push_back( 0.00);
      _eta_bins.push_back( 1.37);
      _eta_bins.push_back( 1.52);
      _eta_bins.push_back( 2.37);

      _eta_bins_areaoffset.push_back(0.0);
      _eta_bins_areaoffset.push_back(1.5);
      _eta_bins_areaoffset.push_back(3.0);
    }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      addProjection(fs, "FS");

      // Consider the final state jets for the energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5);
      /// @todo Memory leak! (Once per run, not serious)
      _area_def = new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()); ///< @todo FastJets should deep copy the pointer if possible
      fj.useJetArea(_area_def);
      addProjection(fj, "KtJetsD05");

      // Consider the leading pt photon with |eta| < 2.37 and pT > 100 GeV
      LeadingParticlesFinalState photonfs(FinalState(-2.37, 2.37, 100.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      addProjection(photonfs, "LeadingPhoton");

      // Book the dsigma/dEt (in eta bins) histograms
      for (size_t i = 0; i < _eta_bins.size() - 1; i++) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue; // skip this bin
        _h_Et_photon[i] = bookHisto1D(1, 1, i+1);
      }

      // Book the dsigma/d|eta| histogram
      _h_eta_photon = bookHisto1D(1,2,1);
    }


    /// Return eta bin for either dsigma/dET histogram (area_eta=false) or energy density correction (area_eta=true)
    size_t _getEtaBin(double eta_w, bool area_eta) const {
      const double eta = fabs(eta_w);
      if (!area_eta) {
        return binIndex(eta, _eta_bins);
      } else {
        return binIndex(eta, _eta_bins_areaoffset);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve leading photon
      Particles photons = applyProjection<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1) vetoEvent;
      const Particle& leadingPhoton = photons[0];

      // Veto events with photon in ECAL crack
      if (inRange(leadingPhoton.abseta(), 1.37, 1.52)) vetoEvent;

      // Compute isolation energy in cone of radius .4 around photon (all particles)
      FourMomentum mom_in_EtCone;
      Particles fs = applyProjection<FinalState>(event, "FS").particles();
      foreach (const Particle& p, fs) {
        // Check if it's outside the cone of 0.4
        if (deltaR(leadingPhoton, p) >= 0.4) continue;
        // Don't count particles in the 5x7 central core
        if (deltaEta(leadingPhoton, p) < .025*5.0*0.5 &&
            deltaPhi(leadingPhoton, p) < (PI/128.)*7.0*0.5) continue;
        // Increment isolation energy
        mom_in_EtCone += p.momentum();
      }

      // Get the area-filtered jet inputs for computing median energy density, etc.
      vector<double> ptDensity, ptSigma, nJets;
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      FastJets fast_jets =applyProjection<FastJets>(event, "KtJetsD05");
      const fastjet::ClusterSequenceArea* clust_seq_area = fast_jets.clusterSeqArea();
      foreach (const Jet& jet, fast_jets.jets()) {
        const double area = clust_seq_area->area(jet);
        if (area > 10e-4 && jet.abseta() < _eta_bins_areaoffset.back())
          ptDensities.at( _getEtaBin(jet.abseta(), true) ).push_back(jet.pT()/area);
      }
      // Compute the median energy density, etc.
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; b++) {
        const int njets = ptDensities[b].size();
        const double ptmedian = (njets > 0) ? median(ptDensities[b]) : 0;
        const double ptsigma = (njets > 0) ? ptDensities[b][(size_t)(0.15865*njets)] : 0;
        nJets.push_back(njets);
        ptDensity.push_back(ptmedian);
        ptSigma.push_back(ptsigma);
      }
      // Compute the isolation energy correction (cone area*energy density)
      const double etCone_area = PI*sqr(0.4) - (7.0*.025)*(5.0*PI/128.);
      const double correction = ptDensity[_getEtaBin(leadingPhoton.abseta(), true)]*etCone_area;

      // Apply isolation cut on area-corrected value
      if (mom_in_EtCone.Et() - correction > 7*GeV) vetoEvent;

      // Fill histograms
      const size_t eta_bin = _getEtaBin(leadingPhoton.abseta(), false);
      _h_Et_photon[eta_bin]->fill(leadingPhoton.Et(), event.weight());
      _h_eta_photon->fill(leadingPhoton.abseta(), event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i = 0; i < _eta_bins.size()-1; i++) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue;
        scale(_h_Et_photon[i], crossSection()/picobarn/sumOfWeights());
      }
      scale(_h_eta_photon, crossSection()/picobarn/sumOfWeights());
    }


  private:

    Histo1DPtr _h_Et_photon[3];
    Histo1DPtr _h_eta_photon;

    fastjet::AreaDefinition* _area_def;

    vector<double> _eta_bins, _eta_bins_areaoffset;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1263495);

}
