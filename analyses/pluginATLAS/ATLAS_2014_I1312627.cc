// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Measurement of V+jets distributions, taken as ratios between W and Z channels
  class ATLAS_2014_I1312627 : public Analysis {
  public:

    /// @name Plotting helper types
    //@{

    struct Plots {
      string ref;
      Histo1DPtr comp[2]; // (de)nominator components
      Scatter2DPtr ratio; // Rjets plot
    };

    //@}


    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2014_I1312627);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // get option
      _mode = 0;
      if ( getOption("LMODE") == "EL" )  _mode = 1;
      if ( getOption("LMODE") == "MU" )  _mode = 2;

      // Set up cuts
      Cut cuts;
      if (_mode == 2) { // muon channel
        cuts = Cuts::pT > 25*GeV && Cuts::abseta < 2.4;;
      } else if (_mode) { // electron channel
        cuts = Cuts::pT > 25*GeV && ( Cuts::etaIn(-2.47, -1.52) || Cuts::etaIn(-1.37, 1.37) || Cuts::etaIn(1.52, 2.47) );
      } else { // combined data extrapolated to common phase space
        cuts = Cuts::pT > 25*GeV && Cuts::abseta < 2.5;
      }

      // Boson finders
      FinalState fs;
      WFinder wfinder(fs, cuts, _mode > 1? PID::MUON : PID::ELECTRON, 40*GeV, 8*TeV, 0., 0.1, 
				              WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, 
                      WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      declare(wfinder, "WF");

      ZFinder zfinder(fs, cuts, _mode > 1? PID::MUON : PID::ELECTRON, 66*GeV, 116*GeV, 0.1, 
                      ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::NO);
      declare(zfinder, "ZF");

      // Jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("WF"));
      jet_fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZF"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL);
      declare(jets, "Jets");


      // Book Rjets plots
      _suff = string("-y0") + to_str(_mode + 1);
      hInit("Njets_incl",  "d01"); // inclusive number of jets
      hInit("Njets_excl",  "d04"); // exclusive number of jets
      hInit("Pt1_N1incl",  "d05"); // leading jet pT, at least 1 jet
      hInit("Pt1_N1excl",  "d06"); // leading jet pT, exactly 1 jet
      hInit("Pt1_N2incl",  "d07"); // leading jet pT, at least 2 jets
      hInit("Pt1_N3incl",  "d08"); // leading jet pT, at least 3 jets
      hInit("Pt2_N2incl",  "d09"); // subleading jet pT, at least 2 jets
      hInit("Pt3_N3incl",  "d10"); // sub-subleading jet pT, at least 3 jets
      hInit("ST_N2incl",   "d11"); // scalar jet pT sum, at least 2 jets
      hInit("ST_N2excl",   "d12"); // scalar jet pT sum, exactly 2 jets
      hInit("ST_N3incl",   "d13"); // scalar jet pT sum, at least 3 jets
      hInit("ST_N3excl",   "d14"); // scalar jet pT sum, exactly 3 jets
      hInit("DR_N2incl",   "d15"); // deltaR(j1, j2), at least 2 jets
      hInit("DPhi_N2incl", "d16"); // deltaPhi(j1, j2), at least 2 jets
      hInit("Mjj_N2incl",  "d17"); // mjj, at least 2 jets
      hInit("Rap1_N1incl", "d18"); // leading jet rapidity, at least 1 jet
      hInit("Rap2_N2incl", "d19"); // subleading jet rapidity, at least 2 jets
      hInit("Rap3_N3incl", "d20"); // sub-subleading jet rapidity, at least 3 jets

      // Also book numerator and denominator for Rjets plots
      for (auto& item : _plots) {
        book(item.second.comp[0], item.second.ref + "2" + _suff, refData(item.second.ref + "1" + _suff));
        book(item.second.comp[1], item.second.ref + "3" + _suff, refData(item.second.ref + "1" + _suff));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve boson candidate
      const WFinder& wf = apply<WFinder>(event, "WF");
      const ZFinder& zf = apply<ZFinder>(event, "ZF");
      if (wf.empty() && zf.empty())  vetoEvent;

      // Retrieve jets
      const JetAlg& jetfs = apply<JetAlg>(event, "Jets");
      Jets jets = jetfs.jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);

      // Apply boson cuts and fill histograms
      if (!zf.empty()) {
        const Particles& leptons = zf.constituents();
        if (oppSign(leptons[0], leptons[1]) && deltaR(leptons[0], leptons[1]) > 0.2)
          fillPlots(leptons, jets, 1);
      }
      if (!wf.empty()) {
        const Particles& leptons = wf.constituentLeptons();
        if (wf.constituentNeutrino().pT() > 25*GeV && wf.mT() > 40*GeV )
          fillPlots(leptons, jets, 0);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      ///  Normalise, scale and otherwise manipulate histograms here
      const double sf( crossSection() / sumOfWeights() );
      for (const auto& item : _plots) {
        scale(item.second.comp[0], sf);
        scale(item.second.comp[1], sf);
        divide(item.second.comp[0], item.second.comp[1], item.second.ratio);
      }
    }
    //@}


    /// Analysis helper functions
    //@{

    /// Histogram filling function, to avoid duplication
    void fillPlots(const Particles& leptons, Jets& jets, int isZ) {
      // Do jet-lepton overlap removal
      idiscardIfAnyDeltaRLess(jets, leptons, 0.5);

      // Calculate jet ST
      const size_t njets = jets.size();
      const double ST = sum(jets, pT, 0.0)/GeV;

      // Fill jet histos
      _plots["Njets_excl"].comp[isZ]->fill(njets + 0.5);
      for (size_t i = 0; i <= njets; ++i)
        _plots["Njets_incl"].comp[isZ]->fill(i + 0.5);

      if (njets > 0) {
        const double pT1  = jets[0].pT()/GeV;
        const double rap1 = jets[0].absrap();
        _plots["Pt1_N1incl" ].comp[isZ]->fill(pT1);
        _plots["Rap1_N1incl"].comp[isZ]->fill(rap1);

        if (njets == 1) {
          _plots["Pt1_N1excl"].comp[isZ]->fill(pT1);
        } else if (njets > 1) {
          const double pT2  = jets[1].pT()/GeV;
          const double rap2 = jets[1].absrap();
          const double dR   = deltaR(jets[0], jets[1]);
          const double dPhi = deltaPhi(jets[0], jets[1]);
          const double mjj  = (jets[0].momentum() + jets[1].momentum()).mass()/GeV;
          _plots["Pt1_N2incl" ].comp[isZ]->fill(pT1);
          _plots["Pt2_N2incl" ].comp[isZ]->fill(pT2);
          _plots["Rap2_N2incl"].comp[isZ]->fill(rap2);
          _plots["DR_N2incl"  ].comp[isZ]->fill(dR);
          _plots["DPhi_N2incl"].comp[isZ]->fill(dPhi);
          _plots["Mjj_N2incl" ].comp[isZ]->fill(mjj);
          _plots["ST_N2incl"  ].comp[isZ]->fill(ST);

          if (njets == 2) {
            _plots["ST_N2excl"].comp[isZ]->fill(ST);
          } else if (njets > 2) {
            const double pT3  = jets[2].pT()/GeV;
            const double rap3 = jets[2].absrap();
            _plots["Pt1_N3incl" ].comp[isZ]->fill(pT1);
            _plots["Pt3_N3incl" ].comp[isZ]->fill(pT3);
            _plots["Rap3_N3incl"].comp[isZ]->fill(rap3);
            _plots["ST_N3incl"  ].comp[isZ]->fill(ST);

            if (njets == 3)
              _plots["ST_N3excl"].comp[isZ]->fill(ST);
          }
        }
      }
    }


    /// Helper for histogram initialisation
    void hInit(string label, string ident) {
      string pre = ident + "-x0";
      _plots[label].ref = pre;
      book(_plots[label].ratio, pre + "1" + _suff, true);
    }

    //@}


  protected:

    // Data members
    size_t _mode;
    string _suff;


  private:

    /// @name Map of histograms
    map<string, Plots> _plots;

  };

  // Hooks for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2014_I1312627);
}
