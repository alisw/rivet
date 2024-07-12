// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief measurement of on-shell ZZ at 13 TeV
  class ATLAS_2017_I1625109 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2017_I1625109);

    /// @name Analysis methods
    //@{

    struct Dilepton {
      Dilepton() {};
      Dilepton(const ParticlePair & _leptons) : leptons(_leptons) {}

      FourMomentum momentum() const {
        return leptons.first.mom() + leptons.second.mom();
      }

      ParticlePair leptons;
    };


    struct Quadruplet {

      vector<DressedLepton> getLeptonsSortedByPt() const {
        vector<DressedLepton> out = { leadingDilepton.leptons.first, leadingDilepton.leptons.second,
                                      subleadingDilepton.leptons.first, subleadingDilepton.leptons.second };
        std::sort(out.begin(), out.end(), cmpMomByPt);
        return out;
      }

      Quadruplet(const Dilepton& dilepton1, const Dilepton& dilepton2) {
        if (dilepton1.momentum().pt() > dilepton2.momentum().pt()) {
          leadingDilepton = dilepton1;
          subleadingDilepton = dilepton2;
        }
        else {
          leadingDilepton = dilepton2;
          subleadingDilepton = dilepton1;
        }
        leptonsSortedByPt = getLeptonsSortedByPt();
      }

      FourMomentum momentum() const {
        return leadingDilepton.momentum() + subleadingDilepton.momentum();
      }

      double distanceFromZMass() const {
        return abs(leadingDilepton.momentum().mass() - Z_mass) + abs(subleadingDilepton.momentum().mass() - Z_mass);
      }

      Dilepton leadingDilepton;
      Dilepton subleadingDilepton;
      vector<DressedLepton> leptonsSortedByPt;
    };

    typedef vector<Quadruplet> Quadruplets;

    typedef std::pair<size_t, size_t> IndexPair;


    vector<IndexPair> getOppositeChargePairsIndices(const vector<DressedLepton>& leptons) {
      vector<IndexPair> indices = {};
      if (leptons.size() < 2) return indices;
      for (size_t i = 0; i < leptons.size(); ++i) {
        for (size_t k = i+1; k < leptons.size(); ++k) {
          const auto charge_i = leptons.at(i).charge();
          const auto charge_k = leptons.at(k).charge();
          if (charge_i == -charge_k) {
            indices.push_back(std::make_pair(i, k));
          }
        }
      }
      return indices;
    }

    bool indicesOverlap(IndexPair a, IndexPair b) {
      return (a.first == b.first || a.first == b.second || a.second == b.first || a.second == b.second);
    }


    bool passesHierarchicalPtRequirements(const Quadruplet& quadruplet) {
      const auto& sorted_leptons = quadruplet.leptonsSortedByPt;
      if (sorted_leptons.at(0).pt() < 20*GeV)  return false;
      if (sorted_leptons.at(1).pt() < 15*GeV)  return false;
      if (sorted_leptons.at(2).pt() < 10*GeV)  return false;
      return true;
    }

    bool passesDileptonMinimumMassRequirement(const Quadruplet& quadruplet) {
      const auto& leptons = quadruplet.leptonsSortedByPt;
      for (const Particle& l1 : leptons) {
        for (const Particle& l2 : leptons) {
          if (isSame(l1, l2)) continue;
          if ((l1.pid() + l2.pid() == 0) && ((l1.mom() + l2.mom()).mass() < 5.0*GeV))  return false;
        }
      }
      return true;
    }

    bool passesLeptonDeltaRRequirements(const Quadruplet& quadruplet) {
      const auto& leptons = quadruplet.leptonsSortedByPt;
      for (const Particle& l1 : leptons) {
        for (const Particle& l2 : leptons) {
          if (isSame(l1, l2))  continue;
          // Any lepton flavour:
          if (deltaR(l1.mom(), l2.mom()) < 0.1)  return false;
          // Different lepton flavour:
          if ((l1.abspid() - l2.abspid() != 0) && (deltaR(l1.mom(), l2.mom()) < 0.2))  return false;
        }
      }
      return true;
    }

    Quadruplets formQuadrupletsByChannel(const vector<DressedLepton>& same_flavour_leptons, vector<IndexPair> indices) {
      Quadruplets quadruplets = {};
      for (size_t i = 0; i <  indices.size(); ++i) {
        for (size_t k = i+1; k <  indices.size(); ++k) {
          const auto& pair_i = indices.at(i);
          const auto& pair_k = indices.at(k);
          if (indicesOverlap(pair_i, pair_k))  continue;
          const auto d1 = Dilepton({same_flavour_leptons.at(pair_i.first), same_flavour_leptons.at(pair_i.second)});
          const auto d2 = Dilepton({same_flavour_leptons.at(pair_k.first), same_flavour_leptons.at(pair_k.second)});
          const auto quadruplet = Quadruplet(d1, d2);
          if (passesHierarchicalPtRequirements(quadruplet)) quadruplets.push_back(quadruplet);
        }
      }
      return quadruplets;
    }

    Quadruplets formQuadrupletsByChannel(const vector<DressedLepton>& electrons, vector<IndexPair> e_indices,
                                         const vector<DressedLepton>& muons,     vector<IndexPair> m_indices) {
      Quadruplets quadruplets = {};
      for (const auto& pair_e : e_indices) {
        for (const auto& pair_m : m_indices) {
          const auto d1 = Dilepton({electrons.at(pair_e.first), electrons.at(pair_e.second)});
          const auto d2 = Dilepton({muons.at(pair_m.first), muons.at(pair_m.second)});
          const auto quadruplet = Quadruplet(d1, d2);
          if (passesHierarchicalPtRequirements(quadruplet)) quadruplets.push_back(quadruplet);
        }
      }
      return quadruplets;
    }


    Quadruplets getQuadruplets(const vector<DressedLepton>& electrons, const vector<DressedLepton>& muons) {
      const auto oc_electrons_indices = getOppositeChargePairsIndices(electrons);
      const auto oc_muons_indices = getOppositeChargePairsIndices(muons);

      const auto electron_quadruplets = formQuadrupletsByChannel(electrons, oc_electrons_indices);
      const auto muon_quadruplets = formQuadrupletsByChannel(muons, oc_muons_indices);
      const auto mixed_quadruplets = formQuadrupletsByChannel(electrons, oc_electrons_indices, muons, oc_muons_indices);

      auto quadruplets = electron_quadruplets;
      quadruplets.insert(quadruplets.end(), muon_quadruplets.begin(), muon_quadruplets.end());
      quadruplets.insert(quadruplets.end(), mixed_quadruplets.begin(), mixed_quadruplets.end());

      return quadruplets;
    }


    Quadruplet selectQuadruplet(const Quadruplets& quadruplets) {
      if (quadruplets.empty()) throw std::logic_error("Expect at least one quadruplet! The user should veto events without quadruplets.");
      Quadruplets sortedQuadruplets = quadruplets;
      std::sort(sortedQuadruplets.begin(), sortedQuadruplets.end(), [](const Quadruplet& a, const Quadruplet& b) {
        return a.distanceFromZMass() < b.distanceFromZMass();
      });
      return sortedQuadruplets.at(0);
    }

    /// Book histograms and initialise projections before the run
    void init() {
      const Cut presel = Cuts::abseta < 5 && Cuts::pT > 100*MeV;
      const FinalState fs(presel);

      // Prompt leptons, photons, neutrinos
      // Excluding those from tau decay
      const PromptFinalState photons(presel && Cuts::abspid == PID::PHOTON, false);
      const PromptFinalState bare_elecs(presel && Cuts::abspid == PID::ELECTRON, false);
      const PromptFinalState bare_muons(presel && Cuts::abspid == PID::MUON, false);

      // Baseline lepton and jet declaration
      const Cut lepton_baseline_cuts = Cuts::abseta < 2.7 && Cuts::pT > 5*GeV;
      const DressedLeptons elecs = DressedLeptons(photons, bare_elecs, 0.1, lepton_baseline_cuts);
      const DressedLeptons muons = DressedLeptons(photons, bare_muons, 0.1, lepton_baseline_cuts);
      declare(elecs, "electrons");
      declare(muons, "muons");

      VetoedFinalState jet_input(fs);
      jet_input.addVetoOnThisFinalState(elecs);
      jet_input.addVetoOnThisFinalState(muons);
      declare(FastJets(jet_input, FastJets::ANTIKT, 0.4), "jets");

      // // Book histograms
      book(_h["pT_4l"], 2, 1, 1);
      book(_h["pT_leading_dilepton"], 8, 1, 1);
      book(_h["pT_subleading_dilepton"], 14, 1, 1);
      book(_h["pT_lepton1"], 20, 1, 1);
      book(_h["pT_lepton2"], 26, 1, 1);
      book(_h["pT_lepton3"], 32, 1, 1);
      book(_h["pT_lepton4"], 38, 1, 1);
      book(_h["absy_4l"], 44, 1, 1);
      book(_h["deltay_dileptons"], 50, 1, 1);
      book(_h["deltaphi_dileptons"], 56, 1, 1);
      book(_h["N_jets"], 62, 1, 1);
      book(_h["N_central_jets"], 68, 1, 1);
      book(_h["N_jets60"], 74, 1, 1);
      book(_h["mass_dijet"], 80, 1, 1);
      book(_h["deltay_dijet"], 86, 1, 1);
      book(_h["scalarpTsum_jets"], 92, 1, 1);
      book(_h["abseta_jet1"], 98, 1, 1);
      book(_h["abseta_jet2"], 104, 1, 1);
      book(_h["pT_jet1"], 110, 1, 1);
      book(_h["pT_jet2"], 116, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(Event const & event) {
      const auto& baseline_electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();
      const auto& baseline_muons = apply<DressedLeptons>(event, "muons").dressedLeptons();

      // Form all possible quadruplets passing hierarchical lepton pT cuts
      const auto quadruplets = getQuadruplets(baseline_electrons, baseline_muons);

      if (quadruplets.empty())  vetoEvent;

      // Select the best quadruplet, the one minimising the total distance from the Z pole mass
      auto const quadruplet = selectQuadruplet(quadruplets);

      // Event selection on the best quadruplet
      if (!passesDileptonMinimumMassRequirement(quadruplet)) vetoEvent;
      if (!passesLeptonDeltaRRequirements(quadruplet)) vetoEvent;
      if (!inRange(quadruplet.leadingDilepton.momentum().mass(), 66*GeV, 116*GeV)) vetoEvent;
      if (!inRange(quadruplet.subleadingDilepton.momentum().mass(), 66*GeV, 116*GeV)) vetoEvent;

      // Select jets
      Jets alljets = apply<JetAlg>(event, "jets").jetsByPt(Cuts::pT > 30*GeV);
      for (const DressedLepton& lep : quadruplet.leptonsSortedByPt)
        ifilter_discard(alljets, deltaRLess(lep, 0.4));
      const Jets jets = alljets;
      const Jets centralJets = filterBy(jets, Cuts::abseta < 2.4);
      const Jets pt60Jets = filterBy(jets, Cuts::pT > 60*GeV);

      const auto& leadingDilepton = quadruplet.leadingDilepton.momentum();
      const auto& subleadingDilepton = quadruplet.subleadingDilepton.momentum();
      
      _h["pT_4l"]->fill((leadingDilepton + subleadingDilepton).pt()/GeV);
      _h["pT_leading_dilepton"]->fill(leadingDilepton.pt()/GeV);
      _h["pT_subleading_dilepton"]->fill(subleadingDilepton.pt()/GeV);
      _h["pT_lepton1"]->fill(quadruplet.leptonsSortedByPt.at(0).pt()/GeV);
      _h["pT_lepton2"]->fill(quadruplet.leptonsSortedByPt.at(1).pt()/GeV);
      _h["pT_lepton3"]->fill(quadruplet.leptonsSortedByPt.at(2).pt()/GeV);
      _h["pT_lepton4"]->fill(quadruplet.leptonsSortedByPt.at(3).pt()/GeV);
      _h["absy_4l"]->fill((leadingDilepton + subleadingDilepton).absrapidity());
      _h["deltay_dileptons"]->fill(fabs(leadingDilepton.rapidity() - subleadingDilepton.rapidity()));
      _h["deltaphi_dileptons"]->fill(deltaPhi(leadingDilepton, subleadingDilepton)/pi);
      _h["N_jets"]->fill(jets.size());
      _h["N_central_jets"]->fill(centralJets.size());
      _h["N_jets60"]->fill(pt60Jets.size());
      
      // If at least one jet present
      if (jets.empty())  vetoEvent;
      _h["scalarpTsum_jets"]->fill(sum(jets, pT, 0.)/GeV);
      _h["abseta_jet1"]->fill(jets.front().abseta());
      _h["pT_jet1"]->fill(jets.front().pt()/GeV);
      
      // If at least two jets present
      if (jets.size() < 2)  vetoEvent; 
      _h["mass_dijet"]->fill((jets.at(0).mom() + jets.at(1).mom()).mass()/GeV);
      _h["deltay_dijet"]->fill(fabs(jets.at(0).rapidity() - jets.at(1).rapidity()));
      _h["abseta_jet2"]->fill(jets.at(1).abseta());
      _h["pT_jet2"]->fill(jets.at(1).pt()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalise histograms to cross section
      const double sf = crossSectionPerEvent() / femtobarn;
      scale(_h, sf);
    }
    //@}

  private:
    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    static constexpr double Z_mass = 91.1876;
    //@}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2017_I1625109);
}
