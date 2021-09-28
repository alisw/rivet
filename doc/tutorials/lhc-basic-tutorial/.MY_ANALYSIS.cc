// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class MY_ANALYSIS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MY_ANALYSIS);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      _lmode = 0; // default accepts either channel
      if ( getOption("LMODE") == "EL" )  _lmode = 1;
      if ( getOption("LMODE") == "MU" )  _lmode = 2;

      _invmode = 0; // _invmode deactivated by default
      if ( getOption("INVMODE") == "STRICT" )  _invmode = 1;
      if ( getOption("INVMODE") == "FIDUCIAL" )  _invmode = 2;

      // full detector acceptance
      const FinalState fs_full(Cuts::abseta < 4.9);

      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      declare(bare_leps, "bare_leps");

      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      DressedLeptons dressed_leps(photons, bare_leps, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);
      declare(dressed_leps, "dressed_leptons");

      // jet collection
      FastJets jets(fs_full, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(jets, "jets");

      // treating muons with |eta| beyond 2.5 as invisible
      VetoedFinalState fs_onlyCentralMuons(fs_full, Cuts::abspid == PID::MUON && Cuts::abseta > 2.5);
      declare(MissingMomentum(_invmode == 1? fs_full : fs_onlyCentralMuons), "invisibles");

      // Book histograms
      vector<double> mll_bins = { 66., 74., 78., 82., 84., 86., 88., 89., 90., 91.,
                                  92., 93., 94., 96., 98., 100., 104., 108., 116. };
      if (_invmode) {
        book(_h["pTmiss"], "pTmiss", 10, 0., 100.);
      }
      else {
        //book(_h["mll"], "mass_ll", mll_bins);
        book(_h["mll_bare"],    "mass_ll_bare",    mll_bins);
        book(_h["mll_dressed"], "mass_ll_dressed", mll_bins);
        book(_h["jets_excl"],  "jets_excl",   6, -0.5,  5.5);
        book(_h["bjets_excl"], "bjets_excl",  3, -0.5,  2.5);
        book(_h["HT"],         "HT",          6,  20., 110.);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& bare_leptons = apply<PromptFinalState>(event, "bare_leps").particles(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);

      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "dressed_leptons").dressedLeptons();

      if (_invmode) {
        if (leptons.empty()) {
          const double pTmiss = apply<MissingMomentum>(event, "invisibles").missingPt()/GeV;
          _h["pTmiss"]->fill(pTmiss);
        }
        vetoEvent;
      }

      if (leptons.size() != 2)  vetoEvent;
      if (leptons[0].pid() != -leptons[1].pid())  vetoEvent; // same flavour, opposite charge

      if (_lmode == 1 && leptons[0].abspid() == PID::MUON)  vetoEvent;
      else if (_lmode == 2 && leptons[0].abspid() == PID::ELECTRON)  vetoEvent;

      const double mll = (leptons[0].mom() + leptons[1].mom()).mass()/GeV;
      _h["mll_dressed"]->fill(mll);

      if (bare_leptons.size() == 2) {
        const double mll_bare = (bare_leptons[0].mom() + bare_leptons[1].mom()).mass()/GeV;
        _h["mll_bare"]->fill(mll_bare);
      }

      if (!inRange(mll, 66*GeV, 116*GeV))  vetoEvent;

      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 10*GeV && Cuts::rap < 4.5);
      idiscardIfAnyDeltaRLess(jets, leptons, 0.4);

      _h["jets_excl"]->fill(jets.size());

      const double HT = sum(jets, pT, 0.0)/GeV;
      _h["HT"]->fill(HT);
      
      size_t bTags = count(jets, hasBTag(Cuts::pT > 5*GeV && Cuts::abseta < 2.5));
      _h["bjets_excl"]->fill(bTags);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection() / sumOfWeights();
      scale(_h, sf);

    }

    //@}

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    size_t _lmode, _invmode;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MY_ANALYSIS);
}
