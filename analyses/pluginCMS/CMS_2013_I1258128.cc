#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  /// CMS Z rapidity measurement
  class CMS_2013_I1258128 : public Analysis {
  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2013_I1258128);


    void init() {
      // Full final state
      const FinalState fs(Cuts::abseta < 5);
      declare(fs, "FS");

      // Z finders for electrons and muons
      Cut cuts = Cuts::abseta < 2.1 && Cuts::pT > 20*GeV;
      const ZFinder zfe(fs, cuts, PID::ELECTRON, 76*GeV, 106*GeV);
      const ZFinder zfm(fs, cuts, PID::MUON, 76*GeV, 106*GeV);
      declare(zfe, "ZFE");
      declare(zfm, "ZFM");

      // Try to get the leading photon
      LeadingParticlesFinalState photonfs(FinalState((Cuts::etaIn(-2.5, 2.5) && Cuts::pT >=  40.0*GeV)));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // Jets
      const FastJets jets(fs, FastJets::ANTIKT, 0.5);
      declare(jets, "JETS");

      // Histograms
      book(_hist1YZ      ,1, 1, 1);
      book(_hist1YJet    ,2, 1, 1);
      book(_hist1YSum    ,3, 1, 1);
book(	  _hist1YDif	,4, 1, 1);
book(	  _hist2YPhoton ,5, 1, 1);
book(	  _hist2YJet	,6, 1, 1);
book(	  _hist2YSum	,7, 1, 1);
book(	  _hist2YDif	,8, 1, 1);
    }


    void makeZCut(const Event& event) {
      // Apply the Z finders and veto if no Z found
      const ZFinder& zfe = apply<ZFinder>(event, "ZFE");
      const ZFinder& zfm = apply<ZFinder>(event, "ZFM");
      if (zfe.empty() && zfm.empty()) vetoEvent;

      // Choose the Z candidate
      const Particles& z = (!zfm.empty()) ? zfm.bosons() : zfe.bosons();
      const Particles& clusteredConstituents = (!zfm.empty()) ? zfm.constituents() : zfe.constituents();

      // Insist that the Z is in a high-pT (boosted) regime
      if (z[0].pT() < 40*GeV) return;

      // Build the jets
      const FastJets& jetfs = apply<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      if (jets.empty()) return;

      // Clean the jets against the lepton candidates with a DeltaR cut of 0.5
      vector<const Jet*> cleanedJets;
      for (const Jet& j : jets) {
        bool isolated = true;
        for (const Particle& p : clusteredConstituents) {
          if (deltaR(p, j) < 0.5) {
            isolated = false;
            break;
          }
        }
        if (isolated) cleanedJets.push_back(&j);
      }
      // Require exactly 1 isolated jet
      if (cleanedJets.size() != 1) return;

      // Fill histos
      const double weight = 1.0;
      const double yz = z[0].rapidity();
      const double yjet = cleanedJets[0]->momentum().rapidity();
      _hist1YZ->fill(fabs(yz), weight);
      _hist1YJet->fill(fabs(yjet), weight);
      _hist1YSum->fill(0.5*fabs(yz + yjet), weight);
      _hist1YDif->fill(0.5*fabs(yz - yjet), weight);
    }


    void makePhotonCut(const Event& event) {
        // Get the photon
        const FinalState& photonfs = apply<FinalState>(event, "LeadingPhoton");
        if (photonfs.particles().size() < 1) return;
        const Particle& photon = photonfs.particles().front();
        if (photon.pT() < 40*GeV) return;
        if (fabs(photon.eta()) > 1.4442 ) return;

      // Build the jets
      const FastJets& jetfs = apply<FastJets>(event, "JETS");
      Jets jets = jetfs.jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      if (jets.empty()) return;

      // Clean the jets against the photon candidate with a DeltaR cut of 0.5
      vector<const Jet*> cleanedJets;
      for (const Jet& j : jets)
        if (deltaR(photon, j) > 0.5)
          cleanedJets.push_back(&j);
      // Require exactly 1 jet
      if (cleanedJets.size() != 1) return;

      // Fill histos
      const double weight = 1.0;
      const double ypho = photon.rapidity();
      const double yjet = cleanedJets[0]->momentum().rapidity();
      _hist2YPhoton->fill(fabs(ypho), weight);
      _hist2YJet->fill(fabs(yjet), weight);
      _hist2YSum->fill(0.5*fabs(ypho + yjet), weight);
      _hist2YDif->fill(0.5*fabs(ypho - yjet), weight);
    }


    void analyze(const Event& event) {
      makeZCut(event);
      makePhotonCut(event);
    }


    void finalize() {
      normalizeByContents(_hist1YZ);
      normalizeByContents(_hist1YJet);
      normalizeByContents(_hist1YSum);
      normalizeByContents(_hist1YDif);
      normalizeByContents(_hist2YPhoton);
      normalizeByContents(_hist2YJet);
      normalizeByContents(_hist2YSum);
      normalizeByContents(_hist2YDif);
    }


    // The CMS normalization in this analysis is that the sum over bin contents
    // is equal to 1. This function normalizes to area = area*bin_width.  /
    // @note This is a strange definition... why?
    void normalizeByContents(Histo1DPtr h) {
      normalize(h, h->bin(0).xWidth());
    }


  private:

    Histo1DPtr _hist1YZ, _hist1YJet, _hist1YSum, _hist1YDif;
    Histo1DPtr _hist2YPhoton, _hist2YJet, _hist2YSum, _hist2YDif;

  };


  // Plugin system hook
  DECLARE_RIVET_PLUGIN(CMS_2013_I1258128);

}
