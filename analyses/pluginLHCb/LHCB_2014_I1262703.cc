// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Study of forward Z + jet production at 7 TeV at LHCb
  /// @author W. Barter, A. Bursche, M. Sirendi (Rivet implementation)
  class LHCB_2014_I1262703 : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2014_I1262703);


    /// Initialise histograms and projections
    void init() {

      // Projections
      const Cut mycut = Cuts::eta >= 2.0 && Cuts::eta <= 4.5 && Cuts::pT > 20*GeV;
      ZFinder zfinder(FinalState(), mycut, PID::MUON, 60*GeV, 120*GeV, 0., ZFinder::ClusterPhotons::NONE);
      declare(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.5);
      declare(jetpro, "Jets");

      // Histograms
      book(_h_jet_pT   , 3, 1, 1);
      book(_h_jet_eta20, 4, 1, 1);
      book(_h_jet_eta10, 4, 1, 2);
      book(_h_Z_y20    , 5, 1, 1);
      book(_h_Z_y10    , 5, 1, 2);
      book(_h_Z_pT20   , 6, 1, 1);
      book(_h_Z_pT10   , 6, 1, 2);
      book(_h_dphi20   , 7, 1, 1);
      book(_h_dphi10   , 7, 1, 2);
      book(_h_dy20     , 8, 1, 1);
      book(_h_dy10     , 8, 1, 2);
    }


    /// Do the analysis
    void analyze(const Event & e) {

      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const Particles leptons = zfinder.constituents();

      const Cut jetSelector = Cuts::eta >= 2.0 && Cuts::eta <= 4.5 && Cuts::pT > 10*GeV;
      const Jets jets = apply<FastJets>(e, "Jets").jetsByPt(jetSelector);

      if (jets.empty()) vetoEvent;

      // Clean the jets against the lepton candidates with a deltaR cut of 0.4
      const Jets cleanedJets = filter_discard(jets, [&](const Jet& j) { return any(leptons, deltaRLess(j, 0.4)); });
      // vector<const Jet*> cleanedJets;
      // for (size_t i = 0; i < jets.size(); i++) {
      //   bool isolated = true;
      //   for (size_t j = 0; j < 2; j++) {
      //     if (deltaR(leptons[j], jets[i]) < 0.4) {
      //       isolated = false;
      //       break;
      //     }
      //   }
      //   if (isolated) cleanedJets.push_back(&jets[i]);
      // }

      // Require at least 1 survivor and note if it is above a 20 GeV jet pT threshold
      if (cleanedJets.empty()) vetoEvent;
      const bool above20 = cleanedJets[0].pT() > 20*GeV;
      const double dphi = deltaPhi(zfinder.boson(), cleanedJets[0]);
      const double drap = zfinder.boson().rap() - cleanedJets[0].rap();

      // Fill histograms
      _h_jet_pT->fill(cleanedJets[0].pT()/GeV);
      _h_jet_eta10->fill(cleanedJets[0].eta());
      _h_Z_y10->fill(zfinder.boson().rap());
      _h_Z_pT10->fill(zfinder.boson().pT()/GeV);
      _h_dphi10->fill(dphi);
      _h_dy10->fill(drap);
      if (above20) {
        _h_jet_eta20->fill(cleanedJets[0].eta());
        _h_Z_y20->fill(zfinder.boson().rap());
        _h_Z_pT20->fill(zfinder.boson().pT()/GeV);
        _h_dphi20->fill(dphi);
        _h_dy20->fill(drap);
      }

    }


    /// Finalize
    void finalize() {
      normalize(_h_jet_pT); normalize(_h_jet_eta20); normalize(_h_jet_eta10); 
      normalize(_h_Z_y20);  normalize(_h_Z_y10);     normalize(_h_Z_pT20); 
      normalize(_h_Z_pT10); normalize(_h_dphi20);    normalize(_h_dphi10); 
      normalize(_h_dy20);   normalize(_h_dy10);
    }


    /// Histograms
    Histo1DPtr _h_jet_pT, _h_jet_eta20, _h_jet_eta10, _h_Z_y20, _h_Z_y10, _h_Z_pT20, _h_Z_pT10, _h_dphi20, _h_dphi10, _h_dy20, _h_dy10;

  };


  RIVET_DECLARE_PLUGIN(LHCB_2014_I1262703);

}
