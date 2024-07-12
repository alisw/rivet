// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {

  /// @brief Z pT over a wide mass range
  class CMS_2022_I2079374 : public Analysis {
    public:

      /// Constructor
      RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2022_I2079374);

      /// Book histograms and initialise projections before the run
      void init() {
        // Get options from the new option system
        // default to combined.
        _mode = 2;
        if (getOption("LMODE") == "EL") {
          _mode = 0;
        } else if (getOption("LMODE") == "MU") {
          _mode = 1;
        } else if (getOption("LMODE") == "EMU") {
          _mode = 2;
        } else {
          MSG_INFO("Assuming a mixed electron/muon sample. Set LMODE=EL/MU if your sample contains a single flavor.");
        }

        FinalState fs;
        PromptFinalState pfs(fs);

        PromptFinalState bareMuons(Cuts::abspid == PID::MUON);
        declare(DressedLeptons(pfs, bareMuons, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, true),
                "muons");

        PromptFinalState bareElectrons(Cuts::abspid == PID::ELECTRON);
        declare(DressedLeptons(pfs, bareElectrons, 0.1, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, true),
                "electrons");

        FastJets jets(fs, FastJets::ANTIKT, 0.4);
        declare(jets, "jets");

        // Histograms
        book(_pt_0jet_mass50_76, 1, 1, 1);
        book(_pt_0jet_mass76_106, 3, 1, 1);
        book(_pt_0jet_mass106_170, 5, 1, 1);
        book(_pt_0jet_mass170_350, 7, 1, 1);
        book(_pt_0jet_mass350_1000, 9, 1, 1);

        book(_pt_1jet_mass50_76, 11, 1, 1);
        book(_pt_1jet_mass76_106, 13, 1, 1);
        book(_pt_1jet_mass106_170, 15, 1, 1);
        book(_pt_1jet_mass170_350, 17, 1, 1);

        book(_phistar_mass50_76, 19, 1, 1);
        book(_phistar_mass76_106, 21, 1, 1);
        book(_phistar_mass106_170, 23, 1, 1);
        book(_phistar_mass170_350, 25, 1, 1);
        book(_phistar_mass350_1000, 27, 1, 1);
      }

      /// Z boson finder.
      /// Note: we don't use the standard ZFinder class in order to stick to
      /// the definition of the publication that is simpler than the ZFinder
      /// algorithm
      /// @param leptons pt-ordered of electron or muon collection to use to build
      /// the Z boson
      std::unique_ptr<Particle> zfinder(const std::vector<DressedLepton> &leptons)
      {
        if (leptons.size() < 2) {
          return nullptr;
        }
        // Leading lepton pT cut
        if (leptons[0].pT() < 25.*GeV) {
          return nullptr;
        }
        if (leptons[0].charge() * leptons[1].charge() > 0) {
          return nullptr;
        }

        std::unique_ptr<Particle> candidate(new Particle(PID::ZBOSON, leptons[0].mom() + leptons[1].mom()));
        if (candidate->mass() < 50*GeV || candidate->mass() > 1000*GeV) {
          return nullptr;
        }
        return candidate;
      }

      /// Perform the per-event analysis
      void analyze(const Event& event)
      {
        // Fetch dressed leptons
        auto muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
        auto electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();

        const std::vector<DressedLepton> *dressedLeptons = nullptr;

        //Look for Z->ee
        std::unique_ptr<Particle> z = zfinder(electrons);
        if (z != nullptr) {
          dressedLeptons = &electrons;
          if (_mode == 1) {
            // Z->ee found, but running with LMODE=MU
            vetoEvent;
          }
        } else if (_mode == 0) {
          // LMODE=EL, but no Z->ee found
          vetoEvent;
        } else { //look for Z->mumu
          z = zfinder(muons);
          if (z.get() != nullptr) {
            dressedLeptons = &muons;
          } else { //no Z boson found
            vetoEvent;
          }
        }

        // Cluster jets
        Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 2.4 && Cuts::pT > 30*GeV);
        idiscardIfAnyDeltaRLess(jets, *dressedLeptons, 0.4);

        //------------------------------------------------

        // Compute phi*
        const Particle& lminus = dressedLeptons->at(0).charge() < 0 ? dressedLeptons->at(0) : dressedLeptons->at(1);
        const Particle& lplus  = dressedLeptons->at(0).charge() < 0 ? dressedLeptons->at(1) : dressedLeptons->at(0);
        double phi_acop = M_PI - deltaPhi(lminus, lplus);
        double costhetastar = tanh( 0.5 * (lminus.eta() - lplus.eta()) );
        double sin2thetastar = (costhetastar > 1) ? 0.0 : (1.0 - sqr(costhetastar));
        double phistar = tan(0.5 * phi_acop) * sqrt(sin2thetastar);

        if (z->mass() > 50. && z->mass() < 76.) {
          _pt_0jet_mass50_76->fill(z->pT());
          _phistar_mass50_76->fill(phistar);
        } else if (z->mass() > 76. && z->mass() < 106.) {
          _pt_0jet_mass76_106->fill(z->pT());
          _phistar_mass76_106->fill(phistar);
        } else if (z->mass() > 106. && z->mass() < 170.) {
          _pt_0jet_mass106_170->fill(z->pT());
          _phistar_mass106_170->fill(phistar);
        } else if (z->mass() > 170. && z->mass() < 350.) {
          _pt_0jet_mass170_350->fill(z->pT());
          _phistar_mass170_350->fill(phistar);
        } else if (z->mass() > 350. && z->mass() < 1000.) {
          _pt_0jet_mass350_1000->fill(z->pT());
          _phistar_mass350_1000->fill(phistar);
        }

        //---------------------- at least one jet --------------------------------
        if (jets.size() < 1) {
          return;
        }

        if (z->mass() > 50. && z->mass() < 76.) {
          _pt_1jet_mass50_76->fill(z->pT());
        } else if (z->mass() > 76. && z->mass() < 106.) {
          _pt_1jet_mass76_106->fill(z->pT());
        } else if (z->mass() > 106. && z->mass() < 170.) {
          _pt_1jet_mass106_170->fill(z->pT());
        } else if (z->mass() > 170. && z->mass() < 350.) {
          _pt_1jet_mass170_350->fill(z->pT());
        }
      }

      /// Calculates the ratio between two histograms
      void calculateRatio(int d, const Histo1DPtr &numerator, const Histo1DPtr &denominator) {
        Scatter2DPtr ratio;
        book(ratio, d, 1, 1, true);

        // Should be the case, but better check
        assert(YODA::Histo1D(*ratio).sameBinning(*numerator));

        // The denominator has a finer binning, so rebin it
        auto rebinned_den = denominator->clone();
        rebinned_den.rebin(numerator->xEdges());

        divide(*numerator, rebinned_den, ratio);
      }

      /// Normalise histograms etc., after the run
      void finalize() {
        double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

        if (_mode == 2)  {
          norm /= 2.;
        }

        scale(_pt_0jet_mass50_76, norm);
        scale(_pt_0jet_mass76_106, norm);
        scale(_pt_0jet_mass106_170, norm);
        scale(_pt_0jet_mass170_350, norm);
        scale(_pt_0jet_mass350_1000, norm);

        scale(_pt_1jet_mass50_76, norm);
        scale(_pt_1jet_mass76_106, norm);
        scale(_pt_1jet_mass106_170, norm);
        scale(_pt_1jet_mass170_350, norm);

        scale(_phistar_mass50_76, norm);
        scale(_phistar_mass76_106, norm);
        scale(_phistar_mass106_170, norm);
        scale(_phistar_mass170_350, norm);
        scale(_phistar_mass350_1000, norm);

        // Calculate the ratios
        calculateRatio(29, _pt_0jet_mass50_76, _pt_0jet_mass76_106);
        calculateRatio(31, _pt_0jet_mass106_170, _pt_0jet_mass76_106);
        calculateRatio(33, _pt_0jet_mass170_350, _pt_0jet_mass76_106);
        calculateRatio(35, _pt_0jet_mass350_1000, _pt_0jet_mass76_106);

        calculateRatio(37, _pt_1jet_mass50_76, _pt_1jet_mass76_106);
        calculateRatio(39, _pt_1jet_mass106_170, _pt_1jet_mass76_106);
        calculateRatio(41, _pt_1jet_mass170_350, _pt_1jet_mass76_106);

        calculateRatio(43, _phistar_mass50_76, _phistar_mass76_106);
        calculateRatio(45, _phistar_mass106_170, _phistar_mass76_106);
        calculateRatio(47, _phistar_mass170_350, _phistar_mass76_106);
        calculateRatio(49, _phistar_mass350_1000, _phistar_mass76_106);
      }

    protected:
      size_t _mode;

    private:
      /// @name Histograms
      Histo1DPtr _pt_0jet_mass50_76,
                 _pt_0jet_mass76_106,
                 _pt_0jet_mass106_170,
                 _pt_0jet_mass170_350,
                 _pt_0jet_mass350_1000;
      Histo1DPtr _pt_1jet_mass50_76,
                 _pt_1jet_mass76_106,
                 _pt_1jet_mass106_170,
                 _pt_1jet_mass170_350;
      Histo1DPtr _phistar_mass50_76,
                 _phistar_mass76_106,
                 _phistar_mass106_170,
                 _phistar_mass170_350,
                 _phistar_mass350_1000;
    };

  RIVET_DECLARE_PLUGIN(CMS_2022_I2079374);
}
