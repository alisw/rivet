// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include <unordered_map>

namespace Rivet {


  /// @brief W+D production in pp at 13 TeV
  class ATLAS_2023_I2628732 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2023_I2628732);

    // Book histograms and initialise before the run
    void init() {
      const FinalState fs;

      // lepton kinematic cuts
      Cut cuts(Cuts::pT > 30*GeV && Cuts::abseta < 2.5);

      // Get photons to dress leptons
      FinalState photons(Cuts::abspid == PID::PHOTON);

      // Get dressed leptons
      IdentifiedFinalState lepids(fs, {{PID::ELECTRON, PID::POSITRON, PID::MUON, PID::ANTIMUON}});
      PromptFinalState leptons(lepids, false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, cuts, true);
      declare(dressedleptons, "DressedLeptons");

      // unstable final-state for Ds
      declare(UnstableParticles(), "UFS");

      // Fiducial cross sections vs species:
      // D+ and D* production fractions can be reweighted to the world average values.
      // We can calculate them for each MC sample separately with the `CharmSpecies`
      // histogram from the Rivet routine, e.g.:
      //   - f(D+) = CharmSpecies->GetBinContent(2) / CharmSpecies->Integral(2,8)
      //   - f(D*) = CharmSpecies->GetBinContent(1) / CharmSpecies->Integral(2,8)
      book(_h["CharmSpecies"], "_CharmSpecies", 8, 0, 8);

      // Differential cross sections per bin
      bookPair("lep_minus", "Dplus", "D_pt",         3);
      bookPair("lep_plus",  "Dplus", "D_pt",         4);
      bookPair("lep_minus", "Dplus", "lep_abs_eta",  5);
      bookPair("lep_plus",  "Dplus", "lep_abs_eta",  6);
      bookPair("lep_minus", "Dstar", "D_pt",         7);
      bookPair("lep_plus",  "Dstar", "D_pt",         8);
      bookPair("lep_minus", "Dstar", "lep_abs_eta",  9);
      bookPair("lep_plus",  "Dstar", "lep_abs_eta", 10);
    }

    /// Perform the per-event analysis
    void analyze(const Event &event) {
      // Retrieve the dressed electrons
      const Particles &signal_leptons = apply<DressedLeptons>(event, "DressedLeptons").particlesByPt();
      if (signal_leptons.size() != 1)  vetoEvent;

      const Particle &lepton = signal_leptons[0];
      const std::string lepton_name = _lepton_names.at(lepton.pid());

      // Get the charm hadrons
      const UnstableParticles &ufs = apply<UnstableFinalState>(event, "UFS");
      std::unordered_map<unsigned int, Particles> particles;

      // Loop over particles
      for (const Particle &p : ufs.particles()) {
        const int id = p.abspid();
        const double pt = p.pT() / GeV;
        const double eta = p.abseta();
        if (_charm_hadron_names.count(id) && pt > 8.0 && eta < 2.2) {
          particles[id].push_back(p);
        }
      }

      // Fill histograms
      for (auto &kv : particles) {
        const unsigned int absPdgId = kv.first;
        const std::string hadron_name = _charm_hadron_names.at(absPdgId);

        for (auto &p : kv.second) {
          // Weight: +1 for OS and -1 for SS
          float charm_charge = (absPdgId == 421) ? p.pid() : p.charge();
          double weight = (charm_charge * lepton.charge() < 0) ? +1.0 : -1.0;

          // Fill charm species for production fraction reweighting
          _h["CharmSpecies"]->fill(_charm_species_map.at(absPdgId), weight);

          // Fill only D+ and D* histograms
          if (absPdgId != PID::DPLUS && absPdgId != PID::DSTARPLUS)  continue;

          // pT(D) overflow
          // Last pT(D) bin extends to infinity (150 only for practical purposes)
          double pt = p.pT() / GeV;
          if (pt >= _max_D_pt)  pt = _max_D_pt - 10;

          // Fill histograms
          _h[histo_name(lepton_name, hadron_name, "lep_abs_eta")]->fill(std::abs(lepton.eta()), weight);
          _h[histo_name(lepton_name, hadron_name, "D_pt")]->fill(pt, weight);
          _h[histo_name(lepton_name, hadron_name, "lep_abs_eta") + "_norm"]->fill(std::abs(lepton.eta()), weight);
          _h[histo_name(lepton_name, hadron_name, "D_pt") + "_norm"]->fill(pt, weight);
        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Finalize
    void finalize() {

      scale(_h, crossSectionPerEvent());

      // D+ and D* production fractions
      const double sum = _h["CharmSpecies"]->integral(false);
      const double fDplus = safediv(_h["CharmSpecies"]->bin(1).sumW(), sum);
      const double fDstar = safediv(_h["CharmSpecies"]->bin(0).sumW(), sum);

      // Reweight to values used in the paper:
      // f(D+) = 0.2404
      // f(D*) = 0.2429
      for (const string lepton_name : {"lep_minus", "lep_plus"}) {
        for (const string hadron_name : {"Dplus", "Dstar"}) {
          const double sf = hadron_name == "Dplus"? (0.2404/fDplus) : (0.2429/fDstar);
          scale(_h[histo_name(lepton_name, hadron_name, "lep_abs_eta")], sf);
          scale(_h[histo_name(lepton_name, hadron_name, "D_pt")], sf);
          normalize(_h[histo_name(lepton_name, hadron_name, "lep_abs_eta") + "_norm"], 1, true);
          normalize(_h[histo_name(lepton_name, hadron_name, "D_pt") + "_norm"], 1, true);
        }
      }

      // The cross-sections from this are analysis are not differential
      for (auto& item : _h) {
        if (item.first != "CharmSpecies")  barchart(item.second, _s[item.first]);
      }
    }


  private:

    string histo_name(const string& lepton, const string& hadron, const string& val) {
      return lepton + "_" + hadron + "_" + val;
    }


    void bookPair(const string& lepton, const string& hadron,
                  const string& val, unsigned int d) {
      // absolute
      string label = histo_name(lepton, hadron, val);
      book(_h[label], "_"+label, refData(d, 1, 1));
      book(_s[label], d, 1, 1);
      // normalised
      label += "_norm";
      book(_h[label], "_"+label, refData(d, 1, 2));
      book(_s[label], d, 1, 2);
    }

    // Mappting for lepton names
    const std::unordered_map<int, std::string> _lepton_names = {
      {PID::ELECTRON, "lep_minus"},
      {PID::POSITRON, "lep_plus"},
      {PID::MUON, "lep_minus"},
      {PID::ANTIMUON, "lep_plus"},
    };

    // Mapping between pdg id an charm hadron names
    const std::unordered_map<unsigned int, std::string> _charm_hadron_names = {
      {PID::DPLUS, "Dplus"},
      {PID::DSTARPLUS, "Dstar"},
      {PID::D0, "Dzero"},
      {PID::DSPLUS, "Ds"},
      {PID::LAMBDACPLUS, "LambdaC"},
      {PID::XI0C, "XiCzero"},
      {PID::XICPLUS, "XiCplus"},
      {PID::OMEGA0C, "OmegaC"},
    };

    // Needed to fill the CharmSpecies histograms
    const std::unordered_map<unsigned int, float> _charm_species_map = {
      {PID::DPLUS, 1.5},
      {PID::DSTARPLUS, 0.5},
      {PID::D0, 2.5},
      {PID::DSPLUS, 3.5},
      {PID::LAMBDACPLUS, 4.5},
      {PID::XI0C, 5.5},
      {PID::XICPLUS, 6.5},
      {PID::OMEGA0C, 7.5},
    };

    // Histogram map
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _s;

    // Last pT(D) bin extends to infinity (150 only for practical purposes)
    double _max_D_pt = 150;
  };


  RIVET_DECLARE_PLUGIN(ATLAS_2023_I2628732);

}
