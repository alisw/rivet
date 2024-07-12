// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Isolated triphotons at 8 TeV
  class ATLAS_2017_I1644367 : public Analysis {
  public:

    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2017_I1644367);

    // Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs, "FS");

      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      PromptFinalState photonfs(Cuts::abspid == PID::PHOTON && Cuts::abseta < 2.37 && Cuts::pT > 15*GeV);
      declare(photonfs, "Photon");

      // Histograms
      book(_h["etg1"]    ,  1, 1, 1);
      book(_h["etg2"]    ,  2, 1, 1);
      book(_h["etg3"]    ,  3, 1, 1);
      book(_h["dphig1g2"],  4, 1, 1);
      book(_h["dphig1g3"],  5, 1, 1);
      book(_h["dphig2g3"],  6, 1, 1);
      book(_h["detag1g2"],  7, 1, 1);
      book(_h["detag1g3"],  8, 1, 1);
      book(_h["detag2g3"],  9, 1, 1);
      book(_h["mg1g2"]   , 10, 1, 1);
      book(_h["mg1g3"]   , 11, 1, 1);
      book(_h["mg2g3"]   , 12, 1, 1);
      book(_h["mg1g2g3"] , 13, 1, 1);

    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Require at least 2 photons in final state
      const Particles photons = apply<PromptFinalState>(event, "Photon").particlesByPt(Cuts::abseta < 1.37 || Cuts::abseta > 1.5);
      if (photons.size() < 3) vetoEvent;

      // Get jets, and corresponding jet areas
      vector<vector<double> > ptDensities(ETA_BINS.size()-1);
      FastJets fastjets = apply<FastJets>(event, "KtJetsD05");
      const auto clust_seq_area = fastjets.clusterSeqArea();
      for (const Jet& jet : fastjets.jets()) {
        const double area = clust_seq_area->area(jet);
        if (area < 1e-3) continue;
        const int ieta = binIndex(jet.abseta(), ETA_BINS);
        if (ieta != -1) ptDensities[ieta].push_back(jet.pT()/area);
      }

      // Compute median jet properties over the jets in the event
      vector<double> ptDensity;
      for (size_t b = 0; b < ETA_BINS.size()-1; ++b) {
        double median = 0.0;
        if (ptDensities[b].size() > 0) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          int nDens = ptDensities[b].size();
          median = (nDens % 2 == 0) ? (ptDensities[b][nDens/2]+ptDensities[b][(nDens-2)/2])/2 : ptDensities[b][(nDens-1)/2];
        }

        ptDensity.push_back(median);
      }

      // Loop over photons and fill vector of isolated ones
      Particles isolated_photons;
      for (const Particle& photon : photons) {
        if (!photon.isPrompt()) continue;

        // Remove photons in ECAL crack region
        const double eta_P = photon.eta();
        const double phi_P = photon.phi();

        // Compute isolation via particles within an R=0.4 cone of the photon
        const Particles fs = apply<FinalState>(event, "FS").particles();
        FourMomentum mom_in_EtCone;
        for (const Particle& p : fs) {
          // Reject if not in cone
          if (deltaR(photon.momentum(), p.momentum()) > 0.4)  continue;
          // Reject if in the 5x7 cell central core
          if (fabs(eta_P - p.eta()) < 0.025 * 5 * 0.5 && fabs(phi_P - p.phi()) < PI/128. * 7 * 0.5)  continue;
          // Sum momentum
          mom_in_EtCone += p.momentum();
        }

        // Now figure out the correction (area*density)
        const double EtCone_area = M_PI*sqr(0.4) - (7*.025)*(5*M_PI/128.); // cone area - central core rectangle
        const double correction = ptDensity[binIndex(fabs(eta_P), ETA_BINS)] * EtCone_area;

        // Discard the photon if there is more than 11 GeV of cone activity
        // NOTE: Shouldn't need to subtract photon itself (it's in the central core)
        if (mom_in_EtCone.Et() - correction > 10*GeV)  continue;
        // Add isolated photon to list
        isolated_photons.push_back(photon);
      }///loop over photons

      // Require at least two isolated photons
      if (isolated_photons.size() < 3) vetoEvent;

      // Select leading pT pair
      sortByPt(isolated_photons);
      const FourMomentum y1 = isolated_photons[0];
      const FourMomentum y2 = isolated_photons[1];
      const FourMomentum y3 = isolated_photons[2];

      // Leading photon should have pT > 40 GeV, subleading > 30 GeV
      if (y1.pT() < 27*GeV)  vetoEvent;
      if (y2.pT() < 22*GeV)  vetoEvent;
      if (y3.pT() < 15*GeV)  vetoEvent;

      // Require the two photons to be separated (dR>0.4)
      if (deltaR(y1,y2) < 0.45)  vetoEvent;
      if (deltaR(y1,y3) < 0.45)  vetoEvent;
      if (deltaR(y2,y3) < 0.45)  vetoEvent;


      const FourMomentum yyy = y1 + y2 + y3;
      const FourMomentum y1y2 = y1 + y2;
      const FourMomentum y1y3 = y1 + y3;
      const FourMomentum y2y3 = y2 + y3;

      const double Myyy = yyy.mass() / GeV;

      const double dPhiy1y2 = mapAngle0ToPi(deltaPhi(y1, y2));
      const double dPhiy1y3 = mapAngle0ToPi(deltaPhi(y1, y3));
      const double dPhiy2y3 = mapAngle0ToPi(deltaPhi(y2, y3));

      const double dEtay1y2 = fabs(y1.eta() - y2.eta());
      const double dEtay1y3 = fabs(y1.eta() - y3.eta());
      const double dEtay2y3 = fabs(y2.eta() - y3.eta());

      if(Myyy < 50.) vetoEvent;

      // Fill histograms

      _h["etg1"]->fill(y1.pT() / GeV);
      _h["etg2"]->fill(y2.pT() / GeV);
      _h["etg3"]->fill(y3.pT() / GeV);

      _h["dphig1g2"]->fill(dPhiy1y2);
      _h["dphig1g3"]->fill(dPhiy1y3);
      _h["dphig2g3"]->fill(dPhiy2y3);

      _h["detag1g2"]->fill(dEtay1y2);
      _h["detag1g3"]->fill(dEtay1y3);
      _h["detag2g3"]->fill(dEtay2y3);

      _h["mg1g2"]->fill(y1y2.mass() / GeV);
      _h["mg1g3"]->fill(y1y3.mass() / GeV);
      _h["mg2g3"]->fill(y2y3.mass() / GeV);

      _h["mg1g2g3"]->fill(Myyy);
    }


    // Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / (femtobarn * sumOfWeights());
      for (auto &hist : _h) {  scale(hist.second, sf); }
    }


  private:

    map<string, Histo1DPtr> _h;
    const vector<double> ETA_BINS = { 0.0, 1.5, 3.0 };

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2017_I1644367);

}
