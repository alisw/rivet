// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// kT splittings in Z events at 8 TeV
  class ATLAS_2017_I1589844 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructors
    ATLAS_2017_I1589844(const string name="ATLAS_2017_I1589844",
                        const string ref_data="ATLAS_2017_I1589844") : Analysis(name) {
      setRefDataName(ref_data);
    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options from the new option system
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;

      const FinalState fs;

      const Cut cuts_mu = (Cuts::pT > 25*GeV) && (Cuts::abseta < 2.4);
      const Cut cuts_el = Cuts::pT > 25*GeV && (Cuts::abseta <= 1.37 || (Cuts::abseta >= 1.52 && Cuts::abseta < 2.47));

      IdentifiedFinalState bare_mu(fs);
      bare_mu.acceptIdPair(PID::MUON);
      IdentifiedFinalState bare_el(fs);
      bare_el.acceptIdPair(PID::ELECTRON);
      const DressedLeptons muons(fs, bare_mu, 0.1, cuts_mu, true);
      const DressedLeptons elecs(fs, bare_el, 0.1, cuts_el, true);
      declare(muons, "muons");
      declare(elecs, "elecs");

      const ChargedFinalState cfs(Cuts::abseta < 2.5 && Cuts::pT > 0.4*GeV);
      VetoedFinalState jet_fs(cfs);
      jet_fs.addVetoOnThisFinalState(muons);
      jet_fs.addVetoOnThisFinalState(elecs);
      declare(FastJets(jet_fs, FastJets::KT, 0.4), "Kt04Jets");
      declare(FastJets(jet_fs, FastJets::KT, 1.0), "Kt10Jets");

      VetoedFinalState jet_fs_all(Cuts::abseta < 2.5 && Cuts::pT > 0.4*GeV);
      jet_fs_all.addVetoOnThisFinalState(muons);
      jet_fs_all.addVetoOnThisFinalState(elecs);
      FastJets jetpro04_all(jet_fs_all, FastJets::KT, 0.4);
      jetpro04_all.useInvisibles();
      declare(jetpro04_all, "Kt04Jets_all");
      FastJets jetpro10_all(jet_fs_all, FastJets::KT, 1.0);
      jetpro10_all.useInvisibles();
      declare(jetpro10_all, "Kt10Jets_all");

      // Histograms with data binning
      _ndij = 8;
      for (size_t i = 0; i < _ndij; ++i) {
        if (_mode == 0 || _mode == 1) {
          string label = "el_d" + to_str(i) + "_kT4";
          _h[label] = bookHisto1D(i + 1, 1, 1);
          _h[label + "_all"] = bookHisto1D(i + 1, 1, 5);
          label = "el_d" + to_str(i) + "_kT10";
          _h[label] = bookHisto1D(i + 1, 1, 3);
          _h[label + "_all"] = bookHisto1D(i + 1, 1, 7);
        }
        if (_mode == 0 || _mode == 2) {
          string label = "mu_d" + to_str(i) + "_kT4";
          _h[label] = bookHisto1D(i + 1, 1, 2);
          _h[label + "_all"] = bookHisto1D(i + 1, 1, 6);
          label = "mu_d" + to_str(i) + "_kT10";
          _h[label] = bookHisto1D(i + 1, 1, 4);
          _h[label + "_all"] = bookHisto1D(i + 1, 1, 8);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {

      // Check we have a Z candidate:
      const vector<DressedLepton>& muons = apply<DressedLeptons>(e, "muons").dressedLeptons();
      const vector<DressedLepton>& elecs = apply<DressedLeptons>(e, "elecs").dressedLeptons();

      bool e_ok = elecs.size() == 2 && muons.empty();
      bool m_ok = elecs.empty() && muons.size() == 2;
      if (_mode == 0 && !e_ok && !m_ok) vetoEvent;
      if (_mode == 1 && !e_ok )  vetoEvent;
      if (_mode == 2 && !m_ok )  vetoEvent;

      string lep_type = elecs.size()? "el_" : "mu_";

      const vector<DressedLepton>& leptons = elecs.size()? elecs : muons;
      if (leptons[0].charge()*leptons[1].charge() > 0) vetoEvent;

      const double dilepton_mass = (leptons[0].momentum() + leptons[1].momentum()).mass();
      if (!inRange(dilepton_mass, 71*GeV, 111*GeV)) vetoEvent;

      const double weight = e.weight();

      // Get kT splitting scales (charged particles only)
      const FastJets& jetpro04 = applyProjection<FastJets>(e, "Kt04Jets");
      const shared_ptr<fastjet::ClusterSequence> seq04 = jetpro04.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq04->n_particles()); ++i) {
        const double dij = sqrt(seq04->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = lep_type + "d" + to_str(i) + "_kT4";
        _h[label]->fill(dij, weight);
      }
      const FastJets& jetpro10 = applyProjection<FastJets>(e, "Kt10Jets");
      const shared_ptr<fastjet::ClusterSequence> seq10 = jetpro10.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq10->n_particles()); ++i) {
        const double dij = sqrt(seq10->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = lep_type + "d" + to_str(i) + "_kT10";
        _h[label]->fill(dij, weight);
      }

      // Get kT splitting scales (all particles)
      const FastJets& jetpro04_all = applyProjection<FastJets>(e, "Kt04Jets_all");
      const shared_ptr<fastjet::ClusterSequence> seq04_all = jetpro04_all.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq04_all->n_particles()); ++i) {
        const double dij = sqrt(seq04_all->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = lep_type + "d" + to_str(i) + "_kT4_all";
        _h[label]->fill(dij, weight);
      }
      const FastJets& jetpro10_all = applyProjection<FastJets>(e, "Kt10Jets_all");
      const shared_ptr<fastjet::ClusterSequence> seq10_all = jetpro10_all.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq10_all->n_particles()); ++i) {
        const double dij = sqrt(seq10_all->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = lep_type + "d" + to_str(i) + "_kT10_all";
        _h[label]->fill(dij, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSectionPerEvent();
      for (auto& kv : _h) scale(kv.second, sf);
    }

    //@}


  protected:

    // Data members like post-cuts event weight counters go here
    size_t _mode, _ndij;


  private:

    // Histograms
    map<string, Histo1DPtr> _h;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1589844);
}
