// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief Measurements of W+W- boson pair production in pp collisions at 13 TeV
  class CMS_2020_I1814328 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2020_I1814328);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance
      const FinalState fs(Cuts::abseta < 4.9);
      const FinalState fsjet4p5(Cuts::abseta < 4.5);
      const FinalState fsjet2p5(Cuts::abseta < 2.5);

      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jet4p5fs(fsjet4p5, FastJets::ANTIKT, 0.4);
      declare(jet4p5fs, "jets4p5");
      FastJets jet2p5fs(fsjet2p5, FastJets::ANTIKT, 0.4);
      declare(jet2p5fs, "jets2p5");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      bare_leps.acceptTauDecays(false);

      // Dress the prompt bare leptons with prompt photons within dR < 0.1,
      // and apply some fiducial cuts on the dressed leptons
      Cut lepton_cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, lepton_cuts);
      declare(dressed_leps, "leptons");

      // Missing momentum
      declare(MissingMomentum(fs), "MET");

      // Book histograms
      book(_h_WW_njets_norm , 2, 1, 1);
      book(_h_WW_mll_norm   , 4, 1, 1);
      book(_h_WW_ptlmax_norm, 5, 1, 1);
      book(_h_WW_ptlmin_norm, 6, 1, 1);
      book(_h_WW_dphill_norm, 7, 1, 1);
      book(_h_WW_njet0      , 8, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Apply a missing-momentum cut
      if (apply<MissingMomentum>(event, "MET").missingPt() < 20*GeV) return;

      // Retrieve dressed leptons, sorted by pT
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "leptons").dressedLeptons();

      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      Jets jets25 = apply<FastJets>(event, "jets4p5").jetsByPt(Cuts::pT > 25*GeV);
      Jets jets30 = apply<FastJets>(event, "jets4p5").jetsByPt(Cuts::pT > 30*GeV);
      Jets jets35 = apply<FastJets>(event, "jets4p5").jetsByPt(Cuts::pT > 35*GeV);
      Jets jets45 = apply<FastJets>(event, "jets4p5").jetsByPt(Cuts::pT > 45*GeV);
      Jets jets60 = apply<FastJets>(event, "jets4p5").jetsByPt(Cuts::pT > 60*GeV);
      Jets jetsNj = apply<FastJets>(event, "jets2p5").jetsByPt(Cuts::pT > 30*GeV);

      // Remove all jets within dR < 0.4 of a dressed lepton
      idiscardIfAnyDeltaRLess(jets25, leptons, 0.4);
      idiscardIfAnyDeltaRLess(jets30, leptons, 0.4);
      idiscardIfAnyDeltaRLess(jets35, leptons, 0.4);
      idiscardIfAnyDeltaRLess(jets45, leptons, 0.4);
      idiscardIfAnyDeltaRLess(jets60, leptons, 0.4);
      idiscardIfAnyDeltaRLess(jetsNj, leptons, 0.4);

      if (leptons.size() == 2 && leptons[0].pid() * leptons[1].pid() < 0) {
        FourMomentum dilCand = leptons[0].momentum() + leptons[1].momentum();
        if (dilCand.mass() > 20*GeV && dilCand.pT() > 30*GeV){
          double ptlmax = leptons[0].pT(); double ptlmin = leptons[1].pT();
          if (ptlmax < ptlmin) {
            ptlmax = leptons[1].pT(); ptlmin = leptons[0].pT();
          }

          if (leptons[0].abspid() != leptons[1].abspid()) {
            _h_WW_njets_norm ->fill(min((double)jetsNj.size()+1, 2.999));
            _h_WW_mll_norm   ->fill(min(dilCand.mass()/GeV, 1499.999));
            _h_WW_ptlmax_norm->fill(min(ptlmax/GeV, 399.999));
            _h_WW_ptlmin_norm->fill(min(ptlmin/GeV, 149.999));
            _h_WW_dphill_norm->fill(deltaPhi(leptons[0], leptons[1]));
          }

          if (jets25.size() == 0) _h_WW_njet0->fill(1.0);
          if (jets30.size() == 0) _h_WW_njet0->fill(2.0);
          if (jets35.size() == 0) _h_WW_njet0->fill(3.0);
          if (jets45.size() == 0) _h_WW_njet0->fill(4.0);
          if (jets60.size() == 0) _h_WW_njet0->fill(5.0);

        }
      }

    }


    /// @todo Replace with barchart()
    void normalizeToSum(Histo1DPtr hist) {
      double sum = 0.;
      for (size_t i = 0; i < hist->numBins(); ++i) {
        sum += hist->bin(i).height();
      }
      scale(hist, 1./sum);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double norm = (sumOfWeights() != 0) ? crossSection()/picobarn/sumOfWeights() : 1.0;

      normalizeToSum(_h_WW_njets_norm );
      normalizeToSum(_h_WW_mll_norm   );
      normalizeToSum(_h_WW_ptlmax_norm);
      normalizeToSum(_h_WW_ptlmin_norm);
      normalizeToSum(_h_WW_dphill_norm);

      scale(_h_WW_njet0, norm);

    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_WW_njets_norm;
    Histo1DPtr _h_WW_mll_norm, _h_WW_ptlmax_norm, _h_WW_ptlmin_norm, _h_WW_dphill_norm;
    Histo1DPtr _h_WW_njet0;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(CMS_2020_I1814328);

}
