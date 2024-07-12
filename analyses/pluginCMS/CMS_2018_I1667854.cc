// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  /// @brief Differential cross section of Z boson production in association with jets at 13 TeV
  ///
  /// @note Code copied from CMS_2015_I1410737 and adapted by P. Gras to CMS-SMP-16-015,
  class CMS_2018_I1667854 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2018_I1667854);


    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system; defaults to combined e+mu
      _mode = 2;
      if ( getOption("LMODE") == "EL" ) _mode = 0;
      if ( getOption("LMODE") == "MU" ) _mode = 1;
      if ( getOption("LMODE") == "EMU" ) _mode = 2;

      // Projections
      FinalState fs;
      VisibleFinalState visfs(fs);
      VetoedFinalState fs_notaudecay(fs);
      fs_notaudecay.addDecayProductsVeto(PID::TAU);
      fs_notaudecay.addDecayProductsVeto(-PID::TAU);

      IdentifiedFinalState bareMuons(fs_notaudecay);
      bareMuons.acceptIdPair(PID::MUON);
      declare(DressedLeptons(fs, bareMuons, /*dRmax = */0.1,
                             Cuts::abseta < 2.4 && Cuts::pT > 20*GeV), "muons");

      IdentifiedFinalState bareElectrons(fs_notaudecay);
      bareElectrons.acceptIdPair(PID::ELECTRON);
      declare(DressedLeptons(fs, bareElectrons, /*dRmax =*/ 0.1,
                             Cuts::abseta < 2.4 && Cuts::pT > 20*GeV), "electrons");

      FastJets jets(visfs, FastJets::ANTIKT, 0.4);
      declare(jets, "jets");

      // Histograms
      book(_h_excmult_jets_tot, 1, 1, 1);
      book(_h_incmult_jets_tot, 2, 1, 1);
      book(_h_zpt1, 3, 1, 1);
      book(_h_leading_jet_pt_tot, 4, 1, 1);
      book(_h_second_jet_pt_tot, 5, 1, 1);
      book(_h_third_jet_pt_tot, 6, 1, 1);
      book(_h_leading_jet_y_tot, 7, 1, 1);
      book(_h_second_jet_y_tot, 8, 1, 1);
      book(_h_third_jet_y_tot, 9, 1, 1);
      book(_h_ht1_tot, 10, 1, 1);
      book(_h_ht2_tot, 11, 1, 1);
      book(_h_ht3_tot, 12, 1, 1);
      book(_h_ptbal1, 13, 1, 1);
      book(_h_ptbal2, 14, 1, 1);
      book(_h_ptbal3, 15, 1, 1);
      book(_h_jzb, 16, 1, 1);
      book(_h_jzb_ptHigh, 17, 1, 1);
      book(_h_jzb_ptLow, 18, 1, 1);
    }


    /// @brief Z boson finder
    ///
    /// @note We don't use the standard ZFinder class in order to stick to
    /// the definition of the publication that is simpler than the ZFinder algorithm.
    ///
    /// @param leptons pt-ordered list of electrons or muons from which to build the Z boson
    std::unique_ptr<Particle> zfinder(const Particles& leptons) {
      if (leptons.size() < 2) return 0;
      if (leptons[0].charge()*leptons[1].charge() > 0) return 0;
      std::unique_ptr<Particle> cand(new Particle(PID::ZBOSON, leptons[0].mom() + leptons[1].mom()));
      if (cand->mass() < 71.*GeV || cand->mass() > 111.*GeV) return 0;
      return cand;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get leptons
      const Particles& muons = apply<DressedLeptons>(event, "muons").particlesByPt();
      const Particles& electrons = apply<DressedLeptons>(event, "electrons").particlesByPt();

      // Look for Z->ee
      std::unique_ptr<Particle> z = zfinder(electrons);
      const Particles* dressedLeptons = 0;
      if (z.get() != nullptr) {
        dressedLeptons = &electrons;
        if (_mode == 1)
          vetoEvent;
      } else { // look for Z->mumu
        z = zfinder(muons);
        if (z.get() != nullptr) {
          dressedLeptons = &muons;
          if (_mode == 0)
            vetoEvent;
        } else { // no Z boson found
          vetoEvent;
        }
      }

      // Cluster jets
      const FastJets& fj = apply<FastJets>(event, "jets");
      const Jets& jets = fj.jetsByPt(Cuts::absrap < 2.4 && Cuts::pT > 30*GeV);

      // Remove jets overlapping with any of the two selected leptons
      Jets goodjets = filter_discard(jets, [dressedLeptons](const ParticleBase& j) {
          return deltaR(j, (*dressedLeptons)[0]) < 0.4 ||  deltaR(j, (*dressedLeptons)[1]) < 0.4;
        });

      // Compute jet pt scalar sum, H_T:
      double ht = sum(goodjets, [](const ParticleBase& j){return j.pT();}, 0.);

      // Fill jet number integral histograms
      _h_excmult_jets_tot->fill(goodjets.size());
      /// @todo Could be better computed by toIntegral transform on exclusive histo
      for (size_t iJet = 0; iJet <= goodjets.size(); iJet++ )
        _h_incmult_jets_tot->fill(iJet);

      if (goodjets.size() < 1) return;

      // Hadronic recoil:
      FourMomentum recoil;
      for (const auto& j: goodjets) {
        recoil += j.momentum();
      }

      // Jet-Z balance = |recoil_T| - |pt(Z)|
      double jzb  = recoil.pT() - z->pT();

      // pT balance:
      double ptbal = (recoil + z->momentum()).pT();

      // Fill leading-jet histograms
      _h_zpt1->fill(z->pT());
      const Jet& j1 = goodjets[0];
      _h_leading_jet_pt_tot->fill(j1.pT()/GeV);
      _h_leading_jet_y_tot->fill(j1.absrapidity());
      _h_ht1_tot->fill(ht/GeV);
      _h_jzb->fill(jzb/GeV);
      if (z->pT() > 50*GeV) {
        _h_jzb_ptHigh->fill(jzb/GeV);
      } else {
        _h_jzb_ptLow->fill(jzb/GeV);
      }
      _h_ptbal1->fill(ptbal/GeV);

      // Fill 2nd jet histograms
      if (goodjets.size() < 2) return;
      const Jet& j2 = goodjets[1];
      _h_second_jet_pt_tot->fill(j2.pT()/GeV);
      _h_second_jet_y_tot->fill(j2.absrapidity());
      _h_ht2_tot->fill(ht/GeV);
      _h_ptbal2->fill(ptbal/GeV);

      // Fill 3rd jet histograms
      if (goodjets.size() < 3) return;
      const Jet& j3 = goodjets[2];
      _h_third_jet_pt_tot->fill(j3.pT()/GeV);
      _h_third_jet_y_tot->fill(j3.absrapidity());
      _h_ht3_tot->fill(ht/GeV);
      _h_ptbal3->fill(ptbal/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalisation factor
      double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;
      // When running in combined mode, need to average to get lepton xsec
      if (_mode == 2) norm /= 2.;

      // MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14)
      //          << std::fixed << std::setprecision(3) << crossSection() << " pb");
      // MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14)
      //          << std::fixed << std::setprecision(3) << numEvents() );
      // MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14)
      //          << std::fixed << std::setprecision(3) << sumOfWeights());
      // MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14)
      //          << std::fixed << std::setprecision(6) << norm);

      scale(_h_excmult_jets_tot, norm);
      scale(_h_incmult_jets_tot, norm);
      scale(_h_zpt1, norm);
      scale(_h_leading_jet_pt_tot, norm);
      scale(_h_second_jet_pt_tot, norm);
      scale(_h_third_jet_pt_tot, norm);
      scale(_h_leading_jet_y_tot, norm);
      scale(_h_second_jet_y_tot, norm);
      scale(_h_third_jet_y_tot, norm);
      scale(_h_ht1_tot, norm);
      scale(_h_ht2_tot, norm);
      scale(_h_ht3_tot, norm);
      scale(_h_ptbal1, norm);
      scale(_h_ptbal2, norm);
      scale(_h_ptbal3, norm);
      scale(_h_jzb, norm);
      scale(_h_jzb_ptHigh, norm);
      scale(_h_jzb_ptLow, norm);
    }


  protected:

    size_t _mode;


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_excmult_jets_tot,  _h_incmult_jets_tot;
    Histo1DPtr _h_leading_jet_pt_tot, _h_second_jet_pt_tot;
    Histo1DPtr _h_third_jet_pt_tot, _h_fourth_jet_pt_tot;
    Histo1DPtr _h_leading_jet_y_tot, _h_second_jet_y_tot;
    Histo1DPtr _h_third_jet_y_tot, _h_fourth_jet_y_tot;
    Histo1DPtr _h_ht1_tot, _h_ht2_tot, _h_ht3_tot, _h_ht4_tot;
    Histo1DPtr _h_ptbal1, _h_ptbal2, _h_ptbal3;
    Histo1DPtr _h_jzb, _h_jzb_ptHigh, _h_jzb_ptLow;
    Histo1DPtr _h_zpt1;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(CMS_2018_I1667854);

}
