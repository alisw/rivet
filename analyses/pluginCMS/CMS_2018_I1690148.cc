#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/MatrixDiag.hh"
#include "Rivet/Tools/fjcontrib/Nsubjettiness.hh"
#include "Rivet/Tools/fjcontrib/EnergyCorrelator.hh"

namespace Rivet {


  /// Measurement of jet substructure observables in $t\bar{t}$ events from $pp$ collisions at 13~TeV
  class CMS_2018_I1690148 : public Analysis {
  public:

    enum Reconstruction { CHARGED=0, ALL=1 };
    enum Observable { MULT=0, PTDS=1, GA_LHA=2, GA_WIDTH=3, GA_THRUST=4, ECC=5, ZG=6, ZGDR=7, NSD=8, TAU21=9, TAU32=10, TAU43=11, C1_00=12, C1_02=13, C1_05=14, C1_10=15, C1_20=16, C2_00=17, C2_02=18, C2_05=19, C2_10=20, C2_20=21, C3_00=22, C3_02=23, C3_05=24, C3_10=25, C3_20=26, M2_B1=27, N2_B1=28, N3_B1=29, M2_B2=30, N2_B2=31, N3_B2=32 };
    enum Flavor { INCL=0, BOTTOM=1, QUARK=2, GLUON=3 };


    /// Minimal constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2018_I1690148);


    /// @name Analysis methods
    //@{


    /// Set up projections and book histograms
    void init() {

      // Cuts
      particle_cut = (Cuts::abseta < 5.0 && Cuts::pT >  0.*GeV);
      lepton_cut   = (Cuts::abseta < 2.4 && Cuts::pT > 15.*GeV);
      jet_cut      = (Cuts::abseta < 2.4 && Cuts::pT > 30.*GeV);

      // Generic final state
      FinalState fs(particle_cut);

      // Dressed leptons
      ChargedLeptons charged_leptons(fs);
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      PromptFinalState prompt_leptons(charged_leptons);
      prompt_leptons.acceptMuonDecays(true);
      prompt_leptons.acceptTauDecays(true);

      PromptFinalState prompt_photons(photons);
      prompt_photons.acceptMuonDecays(true);
      prompt_photons.acceptTauDecays(true);

      // NB. useDecayPhotons=true allows for photons with tau ancestor; photons from hadrons are vetoed by the PromptFinalState;
      DressedLeptons dressed_leptons(prompt_photons, prompt_leptons, 0.1, lepton_cut, true);
      addProjection(dressed_leptons, "DressedLeptons");

      // Projection for jets
      VetoedFinalState fsForJets(fs);
      fsForJets.addVetoOnThisFinalState(dressed_leptons);
      addProjection(FastJets(fsForJets, FastJets::ANTIKT, 0.4,
                             JetAlg::ALL_MUONS, JetAlg::NO_INVISIBLES), "Jets");

      // Booking of histograms
      int d = 0;
      for (int r = 0; r < 2; ++r) { // reconstruction (charged, all)
        for (int o = 0; o < 33; ++o) { // observable
          d += 1;
          for (int f = 0; f < 4; ++f) { // flavor
            char buffer [11];
            sprintf(buffer, "d%02d-x01-y%02d", d, f+1);
            _h[r][o][f] = bookHisto1D(buffer);
          }
        }
      }
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // select ttbar -> lepton+jets
      const vector<DressedLepton>& leptons = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      int nsel_leptons = 0;
      for (const DressedLepton& lepton : leptons) {
        if (lepton.pt() > 26.) nsel_leptons += 1; else vetoEvent; // found veto lepton
      }
      if (nsel_leptons != 1) vetoEvent;

      const Jets all_jets = applyProjection<FastJets>(event, "Jets").jetsByPt(jet_cut);
      if (all_jets.size() < 4) vetoEvent;

      // categorize jets
      int nsel_bjets = 0;
      int nsel_wjets = 0;
      Jets jets[4];
      for (const Jet& jet : all_jets) {
        // check for jet-lepton overlap -> do not consider for selection
        if (deltaR(jet, leptons[0]) < 0.4) continue;

        bool overlap = false;
        bool w_jet   = false;
        for (const Jet& jet2 : all_jets) {
          if (jet.momentum() == jet2.momentum()) continue;
          // check for jet-jet overlap -> do not consider for analysis
          if (deltaR(jet, jet2) < 0.8)
            overlap = true;
          // check for W candidate
          if (jet.bTagged() or jet2.bTagged()) continue;
          FourMomentum w_cand = jet.momentum() + jet2.momentum();
          if (abs(w_cand.mass() - 80.4) < 15.) w_jet = true;
        }

        // count jets for event selection
        if (jet.bTagged()) nsel_bjets += 1;
        if (w_jet) nsel_wjets += 1;

        // jets for analysis
        if (jet.abseta() > 2. or overlap) continue;

        if (jet.bTagged()) {
          jets[BOTTOM].push_back(jet);
        } else if (w_jet) {
          jets[QUARK].push_back(jet);
        } else {
          jets[GLUON].push_back(jet);
        }
      }

      if (nsel_bjets != 2) vetoEvent;
      if (nsel_wjets < 2) vetoEvent;

      // substructure analysis
      // no loop over incl jets -> more loc but faster
      for (int f = 1; f < 4; ++f) {
        for (const Jet& jet : jets[f]) {
          // apply cuts on constituents
          vector<PseudoJet> particles[2];
          for (const Particle& p : jet.particles(Cuts::pT > 1.*GeV)) {
            particles[ALL].push_back( PseudoJet(p.px(), p.py(), p.pz(), p.energy()) );
            if (p.charge3() != 0)
              particles[CHARGED].push_back( PseudoJet(p.px(), p.py(), p.pz(), p.energy()) );
          }

          if (particles[CHARGED].size() == 0) continue;

          // recluster with C/A and anti-kt+WTA
          PseudoJet ca_jet[2];
          JetDefinition ca_def(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
          ClusterSequence ca_charged(particles[CHARGED], ca_def);
          ClusterSequence ca_all(particles[ALL], ca_def);
          ca_jet[CHARGED] = ca_charged.exclusive_jets(1)[0];
          ca_jet[ALL] = ca_all.exclusive_jets(1)[0];

          PseudoJet akwta_jet[2];
          JetDefinition akwta_def(fastjet::antikt_algorithm, fastjet::JetDefinition::max_allowable_R, fastjet::RecombinationScheme::WTA_pt_scheme);
          ClusterSequence akwta_charged(particles[CHARGED], akwta_def);
          ClusterSequence akwta_all(particles[ALL], akwta_def);
          akwta_jet[CHARGED] = akwta_charged.exclusive_jets(1)[0];
          akwta_jet[ALL]     = akwta_all.exclusive_jets(1)[0];

          // calculate observables
          for (int r = 0; r < 2; ++r) {
            int mult = akwta_jet[r].constituents().size();
            // generalized angularities
            _h[r][MULT][INCL]->fill(mult, weight);
            _h[r][MULT][f]->fill(mult, weight);
            if (mult > 1) {
              double ptds = getPtDs(akwta_jet[r]);
              double ga_lha = calcGA(0.5, 1., akwta_jet[r]);
              double ga_width = calcGA(1., 1., akwta_jet[r]);
              double ga_thrust = calcGA(2., 1., akwta_jet[r]);
              _h[r][PTDS][INCL]->fill(ptds, weight);
              _h[r][PTDS][f]->fill(ptds, weight);
              _h[r][GA_LHA][INCL]->fill(ga_lha, weight);
              _h[r][GA_LHA][f]->fill(ga_lha, weight);
              _h[r][GA_WIDTH][INCL]->fill(ga_width, weight);
              _h[r][GA_WIDTH][f]->fill(ga_width, weight);
              _h[r][GA_THRUST][INCL]->fill(ga_thrust, weight);
              _h[r][GA_THRUST][f]->fill(ga_thrust, weight);
            }
            // eccentricity
            if (mult > 3) {
              double ecc = getEcc(akwta_jet[r]);
              _h[r][ECC][INCL]->fill(ecc, weight);
              _h[r][ECC][f]->fill(ecc, weight);
            }
            // N-subjettiness
            if (mult > 2) {
              double tau21 = getTau(2, 1, ca_jet[r]);
              _h[r][TAU21][INCL]->fill(tau21, weight);
              _h[r][TAU21][f]->fill(tau21, weight);
            }
            if (mult > 3) {
              double tau32 = getTau(3, 2, ca_jet[r]);
              _h[r][TAU32][INCL]->fill(tau32, weight);
              _h[r][TAU32][f]->fill(tau32, weight);
            }
            if (mult > 4) {
              double tau43 = getTau(4, 3, ca_jet[r]);
              _h[r][TAU43][INCL]->fill(tau43, weight);
              _h[r][TAU43][f]->fill(tau43, weight);
            }
            // soft drop
            if (mult > 1) {
              vector<double> sd_results = getZg(ca_jet[r]);
              if (sd_results[0] > 0.) {
                _h[r][ZG][INCL]->fill(sd_results[0], weight);
                _h[r][ZG][f]->fill(sd_results[0], weight);
                _h[r][ZGDR][INCL]->fill(sd_results[1], weight);
                _h[r][ZGDR][f]->fill(sd_results[1], weight);
              }
            }
            int nsd = getNSD(0.007, -1., ca_jet[r]);
            _h[r][NSD][INCL]->fill(nsd, weight);
            _h[r][NSD][f]->fill(nsd, weight);
            // C-series energy correlation ratios
            if (mult > 1) {
              double cn_00 = getC(1, 0.0, ca_jet[r]);
              double cn_02 = getC(1, 0.2, ca_jet[r]);
              double cn_05 = getC(1, 0.5, ca_jet[r]);
              double cn_10 = getC(1, 1.0, ca_jet[r]);
              double cn_20 = getC(1, 2.0, ca_jet[r]);
              _h[r][C1_00][INCL]->fill(cn_00, weight);
              _h[r][C1_02][INCL]->fill(cn_02, weight);
              _h[r][C1_05][INCL]->fill(cn_05, weight);
              _h[r][C1_10][INCL]->fill(cn_10, weight);
              _h[r][C1_20][INCL]->fill(cn_20, weight);
              _h[r][C1_00][f]->fill(cn_00, weight);
              _h[r][C1_02][f]->fill(cn_02, weight);
              _h[r][C1_05][f]->fill(cn_05, weight);
              _h[r][C1_10][f]->fill(cn_10, weight);
              _h[r][C1_20][f]->fill(cn_20, weight);
            }
            if (mult > 2) {
              double cn_00 = getC(2, 0.0, ca_jet[r]);
              double cn_02 = getC(2, 0.2, ca_jet[r]);
              double cn_05 = getC(2, 0.5, ca_jet[r]);
              double cn_10 = getC(2, 1.0, ca_jet[r]);
              double cn_20 = getC(2, 2.0, ca_jet[r]);
              _h[r][C2_00][INCL]->fill(cn_00, weight);
              _h[r][C2_02][INCL]->fill(cn_02, weight);
              _h[r][C2_05][INCL]->fill(cn_05, weight);
              _h[r][C2_10][INCL]->fill(cn_10, weight);
              _h[r][C2_20][INCL]->fill(cn_20, weight);
              _h[r][C2_00][f]->fill(cn_00, weight);
              _h[r][C2_02][f]->fill(cn_02, weight);
              _h[r][C2_05][f]->fill(cn_05, weight);
              _h[r][C2_10][f]->fill(cn_10, weight);
              _h[r][C2_20][f]->fill(cn_20, weight);
            }
            if (mult > 3) {
              double cn_00 = getC(3, 0.0, ca_jet[r]);
              double cn_02 = getC(3, 0.2, ca_jet[r]);
              double cn_05 = getC(3, 0.5, ca_jet[r]);
              double cn_10 = getC(3, 1.0, ca_jet[r]);
              double cn_20 = getC(3, 2.0, ca_jet[r]);
              _h[r][C3_00][INCL]->fill(cn_00, weight);
              _h[r][C3_02][INCL]->fill(cn_02, weight);
              _h[r][C3_05][INCL]->fill(cn_05, weight);
              _h[r][C3_10][INCL]->fill(cn_10, weight);
              _h[r][C3_20][INCL]->fill(cn_20, weight);
              _h[r][C3_00][f]->fill(cn_00, weight);
              _h[r][C3_02][f]->fill(cn_02, weight);
              _h[r][C3_05][f]->fill(cn_05, weight);
              _h[r][C3_10][f]->fill(cn_10, weight);
              _h[r][C3_20][f]->fill(cn_20, weight);
            }
            // M/N-series energy correlation ratios
            if (mult > 2) {
              double m2_b1 = getM(2, 1., ca_jet[r]);
              double m2_b2 = getM(2, 2., ca_jet[r]);
              double n2_b1 = getN(2, 1., ca_jet[r]);
              double n2_b2 = getN(2, 2., ca_jet[r]);
              _h[r][M2_B1][INCL]->fill(m2_b1, weight);
              _h[r][M2_B2][INCL]->fill(m2_b2, weight);
              _h[r][N2_B1][INCL]->fill(n2_b1, weight);
              _h[r][N2_B2][INCL]->fill(n2_b2, weight);
              _h[r][M2_B1][f]->fill(m2_b1, weight);
              _h[r][M2_B2][f]->fill(m2_b2, weight);
              _h[r][N2_B1][f]->fill(n2_b1, weight);
              _h[r][N2_B2][f]->fill(n2_b2, weight);
            }
            if (mult > 3) {
              double n3_b1 = getN(3, 1., ca_jet[r]);
              double n3_b2 = getN(3, 2., ca_jet[r]);
              _h[r][N3_B1][INCL]->fill(n3_b1, weight);
              _h[r][N3_B2][INCL]->fill(n3_b2, weight);
              _h[r][N3_B1][f]->fill(n3_b1, weight);
              _h[r][N3_B2][f]->fill(n3_b2, weight);
            }
          }
        }
      }
    }


    void finalize() {
      for (int r = 0; r < 2; ++r) { // reconstruction (charged, all)
        for (int o = 0; o < 33; ++o) { // observable
          for (int f = 0; f < 4; ++f) { // flavor
            normalize(_h[r][o][f], 1.0, false);
          }
        }
      }
    }

    //@}


  private:

    double deltaR(PseudoJet j1, PseudoJet j2) {
      double deta = j1.eta() - j2.eta();
      double dphi = j1.delta_phi_to(j2);
      return sqrt(deta*deta + dphi*dphi);
    }


    double getPtDs(PseudoJet jet) {
      double mult   = jet.constituents().size();
      double sumpt  = 0.; // would be jet.pt() in WTA scheme but better keep it generic
      double sumpt2 = 0.;
      for (auto p : jet.constituents()) {
        sumpt  += p.pt();
        sumpt2 += pow(p.pt(), 2);
      }
      double ptd = sumpt2/pow(sumpt,2);
      return max(0., sqrt((ptd-1./mult) * mult/(mult-1.)));
    }


    double calcGA(double beta, double kappa, PseudoJet jet) {
      double sumpt = 0.;
      for (const auto& p : jet.constituents()) {
        sumpt += p.pt();
      }
      double ga = 0.;
      for (auto p : jet.constituents()) {
        ga += pow(p.pt()/sumpt, kappa) * pow(deltaR(jet, p)/0.4, beta);
      }
      return ga;
    }


    double getEcc(PseudoJet jet) {
      // Covariance matrix
      Matrix<2> M;
      for (const auto& p : jet.constituents()) {
        Matrix<2> MPart;
        MPart.set(0, 0, (p.eta() - jet.eta()) * (p.eta() - jet.eta()));
        MPart.set(0, 1, (p.eta() - jet.eta()) * mapAngleMPiToPi(p.phi() - jet.phi()));
        MPart.set(1, 0, mapAngleMPiToPi(p.phi() - jet.phi()) * (p.eta() - jet.eta()));
        MPart.set(1, 1, mapAngleMPiToPi(p.phi() - jet.phi()) * mapAngleMPiToPi(p.phi() - jet.phi()));
        M += MPart * p.e();
      }
      // Calculate eccentricity from eigenvalues
      const EigenSystem<2> eigen = diagonalize(M);
      return 1. - eigen.getEigenValues()[1]/eigen.getEigenValues()[0];
    }


    double getTau(int N, int M, PseudoJet jet) {
      fjcontrib::Nsubjettiness::NsubjettinessRatio tau_ratio(N, M, fjcontrib::Nsubjettiness::OnePass_WTA_KT_Axes(),
                                                             fjcontrib::Nsubjettiness::NormalizedMeasure(1.0, 0.4));
      return tau_ratio(jet);
    }


    vector<double> getZg(PseudoJet jet) {
      PseudoJet jet0 = jet;
      PseudoJet jet1, jet2;
      double zg = 0.;
      while (zg < 0.1 && jet0.has_parents(jet1, jet2)) {
        zg   = jet2.pt()/jet0.pt();
        jet0 = jet1;
      }
      if (zg < 0.1) return {-1., -1.};
      return {zg, jet1.delta_R(jet2)};
    }


    int getNSD(double zcut, double beta, PseudoJet jet) {
      PseudoJet jet0 = jet;
      PseudoJet jet1, jet2;
      int nsd = 0.;
      double zg = 0.;
      while (jet0.has_parents(jet1, jet2)) {
        zg = jet2.pt()/jet0.pt();
        if (zg > zcut * pow(jet1.delta_R(jet2)/0.4, beta))
          nsd += 1;
        jet0 = jet1;
      }
      return nsd;
    }


    double getC(int N, double beta, PseudoJet jet) {
      fjcontrib::EnergyCorrelatorDoubleRatio C(N, beta);
      return C(jet);
    }


    double getM(int N, double beta, PseudoJet jet) {
      fjcontrib::EnergyCorrelatorMseries CM(N, beta);
      return CM(jet);
    }


    double getN(int N, double beta, PseudoJet jet) {
      fjcontrib::EnergyCorrelatorNseries CN(N, beta);
      return CN(jet);
    }


  private:

    // @name Histogram data members
    //@{

    Cut particle_cut, lepton_cut, jet_cut;
    Histo1DPtr _h[2][33][4];

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2018_I1690148);

}
