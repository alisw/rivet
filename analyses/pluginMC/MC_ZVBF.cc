// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for Zjj events
  class MC_ZVBF : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_ZVBF);

    /// @name Analysis methods
    //@{

    /// Initialize
    void init() {
      _dR=0.1;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
      _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      FinalState fs;
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      ZFinder zfinder(fs, cut, _lepton, 65*GeV, 115*GeV, _dR, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      const double sqrts = sqrtS() ? sqrtS() : 14*TeV;
      book(_h["gap_inc"], "N_gapjets_inclusive", 8, -0.5, 7.5);
      book(_h["gap_exc"], "N_gapjets_exclusive", 8, -0.5, 7.5);
      book(_h["Z_jet1_deta"], "Z_jet1_deta", 50, -5, 5);
      book(_h["Z_jet1_dR"], "Z_jet1_dR", 25, 0.5, 7.0);
      book(_h["HT"], "jets_HT", logspace(40, 50, sqrts/GeV/2.0));
      book(_h["mjj"], "m_jj", 40, 200.0, sqrts/GeV/2.0);
      book(_h["jve_mjj"], "_jve_mjj", 40, 200.0, sqrts/GeV/2.0);
      book(_h["pTV"] ,"Z_pT", logspace(100, 1.0, 0.5*sqrts/GeV));
      book(_h["dphi"], "dphi_jj", 20., -1., 1.);
      book(_h["drap"], "drap_jj", 20., -10., 10.);
      book(_h["3JC"], "jet_3_centrality", 25., -2.5, 2.5);

      book(_s["jve_mjj"], "jet_veto_efficiency_mjj");

      for (size_t i = 0; i < 4; ++i) {
        const string pTname = "jet_pT_" + to_str(i+1);
        const double pTmax = 1.0/(double(i)+2.0) * sqrts/GeV/2.0;
        const int nbins_pT = 100/(i+1);
        if (pTmax > 10) { // Protection aginst logspace exception, needed for LEP
          book(_h[pTname], pTname, logspace(nbins_pT, 10.0, pTmax));
        }
        const string etaname = "jet_eta_" + to_str(i+1);
        book(_h[etaname], etaname, (i > 1 ? 25 : 50), -5.0, 5.0);
        const string rapname = "jet_y_" + to_str(i+1);
        book(_h[rapname], rapname, (i > 1 ? 25 : 50), -5.0, 5.0);
        const string phiname = "jet_phi_" + to_str(i+1);
        book(_h[phiname], phiname, (i > 1 ? 25 : 50), -1.0, 1.0);
      }
    }


    /// Do the analysis
    void analyze(const Event & e) {
      MSG_TRACE("MC_ZVBF: running ZFinder");
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const FourMomentum& zmom = zfinder.bosons()[0].momentum();
      MSG_TRACE("MC_ZVBF: have exactly one Z boson candidate");

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(Cuts::absrap < 5 && Cuts::pT > 30*GeV);
      if (jets.size() < 2) {
        MSG_TRACE("MC_ZVBF: does not have at least two valid jets");
        vetoEvent;
      }

      Jet tag1 = jets.at(0);
      Jet tag2 = jets.at(1);
      const double mjj = (tag1.mom() + tag2.mom()).mass()/GeV;
      if (mjj < 200.) {
        MSG_TRACE("MC_ZVBF: should have at least 200 GeV in Mjj");
        vetoEvent;
      }

      // jet kinematics
      for (size_t i = 0; i < min(4u, jets.size()); ++i) {
        const string pTname  = "jet_pT_"  + to_str(i+1);
        const string etaname = "jet_eta_" + to_str(i+1);
        const string rapname = "jet_y_"   + to_str(i+1);
        const string phiname = "jet_phi_" + to_str(i+1);
        _h[pTname]->fill(jets[i].pT()/GeV);
        _h[etaname]->fill(jets[i].eta());
        _h[rapname]->fill(jets[i].rap());
        _h[phiname]->fill(mapAngleMPiToPi(jets[i].phi())/M_PI);
      }

      size_t n_gap = 0;
      // start loop at the 3rd hardest pT jet
      for (size_t i = 2; i < jets.size(); ++i) {
        const Jet j = jets.at(i);
        if (isBetween(j, tag1, tag2))  ++n_gap;
      }

      // gap-jet multiplicities
      _h["gap_exc"]->fill(n_gap);
      for (size_t i = 0; i <= 7; ++i) {
        if (n_gap >= i) {
          _h["gap_inc"]->fill(i);
        }
      }

      _h["jve_mjj"]->fill(mjj);
      if (n_gap) {
        // third-jet centrality
        const double rap1 = jets[0].rap();
        const double rap2 = jets[1].rap();
        const double rap3 = jets[2].rap();
        const double JC = (rap3 - 0.5*(rap1 + rap2))/(rap1 - rap2);
        _h["3JC"]->fill(JC);
      }
      else {
        MSG_TRACE("MC_ZVBF: should satisfy a CJV");
        const double HT = sum(jets, Kin::pT, 0.0)/GeV;
        _h["HT"]->fill(HT);
        _h["mjj"]->fill(mjj);
        _h["pTV"]->fill(zmom.pT()/GeV);
        _h["dphi"]->fill(signedDeltaPhi(tag1, tag2));
        _h["drap"]->fill(tag1.rap() - tag2.rap());

        // Z tagging jet correlations
        _h["Z_jet1_deta"]->fill(zmom.eta()-jets[0].eta());
        _h["Z_jet1_dR"]->fill(deltaR(zmom, jets[0].momentum()));
      }
    }


    /// Finalize
    void finalize() {
      scale(_h, crossSection()/femtobarn/sumOfWeights());
      efficiency(_h["mjj"], _h["jve_mjj"], _s["jve_mjj"]);
    }

    //@}

    // check if jet is between tagging jets
    bool isBetween(const Jet &probe, const Jet &boundary1, const Jet &boundary2) {
      double y_p = probe.rapidity();
      double y_b1 = boundary1.rapidity();
      double y_b2 = boundary2.rapidity();

      double y_min = std::min(y_b1, y_b2);
      double y_max = std::max(y_b1, y_b2);

      return  (y_p > y_min && y_p < y_max);
    }

    double signedDeltaPhi(Jet &j1, Jet &j2) {
      double dphijj = 0.;
      if (j1.rap() > j2.rap())  dphijj = j1.phi() - j2.phi();
      else                      dphijj = j2.phi() - j1.phi();
      return mapAngleMPiToPi(dphijj)/M_PI;
    }

  private:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    //@{
    double _dR;
    PdgId _lepton;
    //@}

    /// @name Histograms
    //@{
    map<string,Histo1DPtr> _h;
    map<string,Scatter2DPtr> _s;
    //@}

  };

  RIVET_DECLARE_PLUGIN(MC_ZVBF);
}
