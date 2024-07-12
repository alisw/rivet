// -*- C++ -*-
#ifndef RIVET_RIVETHTT_HH
#define RIVET_RIVETHTT_HH

#include "Rivet/Jet.hh"
#include "HEPTopTagger/HEPTopTagger.hh"

namespace Rivet {


  class HTT {
  public:

    /// Clustering algorithms
    ///
    /// @todo Merge with JetAlg alg enum
    enum Algo { KT=0, AKT=1, ANTIKT=1, CA=2, CAMBRIDGE=2, CAMBRIDGE_AACHEN=2 };


    /// HTT operating mode
    enum Mode {
      EARLY_MASSRATIO_SORT_MASS,     // applies 2D mass plane requirements then select the candidate which minimizes |m_cand-mt|
      LATE_MASSRATIO_SORT_MASS,      // selects the candidate which minimizes |m_cand-mt|
      EARLY_MASSRATIO_SORT_MODDJADE, // applies the 2D mass plane requirements then select the candidate with highest jade distance
      LATE_MASSRATIO_SORT_MODDJADE,  // selects the candidate with highest modified jade distance
      TWO_STEP_FILTER                // only analyzes the candidate built with the highest pT(t) after unclustering
    };


    struct InputParameters {

      /// HTT execution mode
      Mode mode = Mode::EARLY_MASSRATIO_SORT_MASS;

      /// @name Optimal-R parameters
      /// @{
      /// Use the optimal-R finder (default), or use fixed R
      bool do_optimalR = true;
      double optimalR_min = 0.5; // min jet size
      double optimalR_step = 0.1; // step size
      double optimalR_threshold = 0.2; // step size
      /// @}

      /// @name Mass-drop declustering
      /// @{
      double mass_drop = 0.8;
      double max_subjet_mass = 30*GeV;
      /// @}

      /// @name Filtering
      /// @{
      unsigned int filt_N = 5; // set_nfilt
      double filtering_R = 0.3; // max subjet distance for filtering
      double filtering_minpt = 0.; // min subjet pt for filtering
      /// @}

      /// Jet algorithm used for filtering
      Algo filtering_algorithm = Algo::CA;

      /// Reclustering jet algorithm
      Algo reclustering_algorithm = Algo::CA;

      /// @name Top-mass ranges
      /// @{
      /// @todo Take from a central set of constants
      double top_mass = 172.3*GeV;
      /// @todo Take from a central set of constants
      double W_mass = 80.4*GeV;
      double Mtop_min = 150*GeV;
      double Mtop_max = 200*GeV; //set_top_range(min,max)
      /// @}

      /// @name Top-mass ratio range
      /// @{
      double fw = 0.15;
      double mass_ratio_range_min = (1.-fw)*W_mass/top_mass;
      double mass_ratio_range_max = (1.+fw)*W_mass/top_mass;
      /// @}

      /// @name Mass-ratio cuts
      /// @{
      double m23cut = 0.35;
      double m13cutmin = 0.2;
      double m13cutmax = 1.3;
      /// @}

      /// @name Pruning cuts
      /// @{
      double prune_zcut = 0.1;
      double prune_rcut = 0.5;
      /// @}

    };



    /// Constructor without arguments
    HTT() {}

    /// Constructor with arguments
    HTT(HTT::InputParameters& params) {
      setParams(params);
    }

    // /// Destructor
    // ~HTT() {}

    /// Set the tagging parameters
    void setParams(HTT::InputParameters& params);

    /// Run the top tagger on a given jet
    void calc(Jet& jet);

    /// Top jet
    const Jet topJet() const;
    /// The bottom jet inside the top
    const Jet bJet() const;
    /// The W jet inside the top
    const Jet wJet() const;
    /// Leading subjet from W
    const Jet w1Jet() const;
    /// Second leading subjet from W
    const Jet w2Jet() const;

    /// Top jet, as a pseudojet
    const PseudoJet& topJet() const;
    /// The bottom jet inside the top, as a pseudojet
    const PseudoJet& bJet() const;
    /// The W jet inside the top, as a pseudojet
    const PseudoJet& wJet() const;
    /// Leading subjet from W, as a pseudojet
    const PseudoJet& w1Jet() const;
    /// Second leading subjet from W, as a pseudojet
    const PseudoJet& w2Jet() const;

    /// pT-ordered subjets
    const Jets& subjets() const;

    // /// Print tagger information
    // void info() const;
    // /// Print tagger settings
    // void settings() const;

    /// The pruned mass
    double prunedMass() const;

    // The unfiltered mass
    double unfilteredMass() const;

    /// Difference between the reco top mass and the true top mass
    double deltaTopMass() const;

    /// Is the jet tagged?
    bool isTopTagged() const;

    /// Was the top-mass window requirement passed?
    bool passedMassCutTop() const;

    /// 2D mass plane requirements passed?
    bool passedMassCut2D() const;


  private:

    fastjet::HEPTopTagger::HEPTopTagger _tagger;

  };


  /// Below can be moved to HTT.cc at some point

  void HTT::setParams(HTT::InputParameters& params) {
    _tagger = fastjet::HEPTopTagger::HEPTopTagger();

    // Optimal R
    _tagger.do_optimalR(params.do_optimalR);
    _tagger.set_optimalR_min(params.optimalR_min);
    _tagger.set_optimalR_step(params.optimalR_step);
    _tagger.set_optimalR_threshold(params.optimalR_threshold);

    // Candidate selection
    fastjet::HEPTopTagger::Mode mode;
    if (params.mode == HTT::EARLY_MASSRATIO_SORT_MASS) {
      mode = fastjet::HEPTopTagger::EARLY_MASSRATIO_SORT_MASS;
    } else if (params.mode == HTT::LATE_MASSRATIO_SORT_MASS) {
      mode = fastjet::HEPTopTagger::LATE_MASSRATIO_SORT_MASS;
    } else if (params.mode == HTT::EARLY_MASSRATIO_SORT_MODDJADE) {
      mode = fastjet::HEPTopTagger::EARLY_MASSRATIO_SORT_MODDJADE;
    } else if (params.mode == HTT::LATE_MASSRATIO_SORT_MODDJADE) {
      mode = fastjet::HEPTopTagger::LATE_MASSRATIO_SORT_MODDJADE;
    } else {
      mode = fastjet::HEPTopTagger::TWO_STEP_FILTER;
    }
    _tagger.set_mode(mode);
    _tagger.set_mt(params.top_mass);
    _tagger.set_mw(params.W_mass);
    _tagger.set_top_mass_range(params.Mtop_min, params.Mtop_max);
    _tagger.set_fw(params.fw);
    _tagger.set_mass_ratio_range(params.mass_ratio_range_min, params.mass_ratio_range_max);
    _tagger.set_mass_ratio_cut(params.m23cut, params.m13cutmin, params.m13cutmax);

    // Filtering
    _tagger.set_filtering_n(params.filt_N);
    _tagger.set_filtering_R(params.filtering_R);
    _tagger.set_filtering_minpt_subjet(params.filtering_minpt);

    fastjet::JetAlgorithm algo = fastjet::antikt_algorithm;
    if (params.filtering_algorithm == Algo::CA) algo = fastjet::cambridge_algorithm;
    else if (params.filtering_algorithm == Algo::KT) algo = fastjet::kt_algorithm;
    _tagger.set_filtering_jetalgorithm(algo);

    // Reclustering
    algo = fastjet::antikt_algorithm;
    if (params.reclustering_algorithm == Algo::CA) algo = fastjet::cambridge_algorithm;
    else if (params.reclustering_algorithm == Algo::KT) algo = fastjet::kt_algorithm;
    _tagger.set_reclustering_jetalgorithm(algo);

    // Mass-drop
    _tagger.set_mass_drop_threshold(params.mass_drop);
    _tagger.set_mass_drop_threshold(params.mass_drop);

    // Pruning
    _tagger.set_pruning_rcut_factor(params.prune_rcut);
    _tagger.set_pruning_zcut(params.prune_zcut);
  }


  void HTT::calc(Jet& jet) {
    _tagger.run(jet.pseudojet());
  }


  const Jet HTT::topJet() const { return Jet(topPjet()); }
  const Jet HTT::bJet() const { return Jet(bPjet()); }
  const Jet HTT::wJet() const { return Jet(wPjet()); }
  const Jet HTT::w1Jet() const { return Jet(w1Pjet()); }
  const Jet HTT::w2Jet() const { return Jet(w2Pjet()); }


  const PseudoJet& HTT::topPjet() const { return _tagger.t(); }
  const PseudoJet& HTT::bPjet() const { return _tagger.b(); }
  const PseudoJet& HTT::wPjet() const { return _tagger.W(); }
  const PseudoJet& HTT::w1Pjet() const { return _tagger.W1(); }
  const PseudoJet& HTT::w2Pjet() const { return _tagger.W2(); }


  Jets HTT::subjets() const {
    Jets rtn;
    rtn.reserve(3);
    rtn.emplace_back(_tagger.j1());
    rtn.emplace_back(_tagger.j2());
    rtn.emplace_back(_tagger.j3());
    return rtn;
  }


  // void HTT::info() const { _tagger.get_info(); }
  // void HTT::settings() const { _tagger.get_setting(); }

  double HTT::prunedMass() const { return _tagger.pruned_mass(); }

  double HTT::unfilteredMass() const { return _tagger.unfiltered_mass(); }

  double HTT::deltaTopMass() const { return _tagger.delta_top(); }

  bool HTT::isTopTagged() const { return _tagger.is_tagged(); }

  bool HTT::passedMassCutTop() const { return _tagger.is_maybe_top(); }

  bool HTT::passedMassCut2D() const { return _tagger.is_masscut_passed(); }

  /// Above can be moved to RivetHTT.cc at some point


}

#endif
