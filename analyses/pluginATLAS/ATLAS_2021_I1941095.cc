// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/MendelMin.hh"

#include "fastjet/tools/Filter.hh"

namespace Rivet {


  /// @brief Energy asymmetry in ttj at 13 TeV
  class ATLAS_2021_I1941095 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2021_I1941095);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Declare projections

      // Photons
      PromptFinalState promptphotons(Cuts::abspid == PID::PHOTON, false);

      // Electrons
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true); // true = use electrons from prompt tau decays
      DressedLeptons all_dressed_el(promptphotons, bare_el, 0.1, Cuts::abseta < 2.5, true);
      DressedLeptons electrons(promptphotons, bare_el, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25*GeV, true);
      declare(electrons,"electrons");

      // Muons
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true); // true = use muons from prompt tau decays
      DressedLeptons all_dressed_mu(promptphotons, bare_mu, 0.1, Cuts::abseta < 2.5, true);
      DressedLeptons muons(promptphotons,bare_mu, 0.1, Cuts::abseta <2.5 && Cuts::pT > 25*GeV, true);
      declare(muons,"muons");

      // AntiKt4TruthWZJets as AntiKt4TruthWZJets, but w/o photons from hadrons in dressing
      const InvisibleFinalState invisibles(true, true);
      VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0)); // changed from 4.5 to 5.0
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);
      vfs.addVetoOnThisFinalState(invisibles); // new
      FastJets jets(vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::ALL); // changed invisible from DECAY to ALL
      declare(jets,"jets");

      // AntiKt10TruthTrimmedPtFrac5SmallR20Jets
      FinalState fs(Cuts::abseta < 5.0);
      FastJets fjets(fs, FastJets::ANTIKT, 1.0, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      _trimmer = fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));
      declare(fjets,"fjets");

      // Missing momentum
      declare(MissingMomentum(), "MissingMomentum");

      // Parton level top quarks after FSR
      // options are: decaymode, emu_from_prompt_tau, include_hadronic_taus
      declare(PartonicTops(PartonicTops::DecayMode::E_MU, true, false), "PartonicTops_EMU");
      declare(PartonicTops(PartonicTops::DecayMode::E_MU, false, false), "PartonicTops_EMU_notau");
      declare(PartonicTops(PartonicTops::DecayMode::HADRONIC, false, true), "PartonicTops_HADRONIC");
      declare(PartonicTops(PartonicTops::DecayMode::HADRONIC, false, false), "PartonicTops_HADRONIC_notau");

      // Book histograms
      const Scatter2D& ref_asymm = refData(1, 1, 1);
      book(_h["pos"], "_thetaj_opt_depos", ref_asymm);
      book(_h["neg"], "_thetaj_opt_deneg", ref_asymm);
      book(_asymm, 1, 1, 1, true);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Parton-level top quarks // after FSR
      const Particles partonicTops_EMU = apply<ParticleFinder>(event, "PartonicTops_EMU").particlesByPt();
      const Particles partonicTops_EMU_notau = apply<ParticleFinder>(event, "PartonicTops_EMU_notau").particlesByPt();
      const Particles partonicTops_HADRONIC = apply<ParticleFinder>(event, "PartonicTops_HADRONIC").particlesByPt();
      const Particles partonicTops_HADRONIC_notau = apply<ParticleFinder>(event, "PartonicTops_HADRONIC_notau").particlesByPt();

      // Filter semi-leptonic (e,mu,tau) events: Veto dileptonic/nonleptonic events and events with 2 taus
      int nLeptons = partonicTops_EMU.size() + partonicTops_HADRONIC.size() - partonicTops_HADRONIC_notau.size();
      if (nLeptons != 1) vetoEvent;

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();
      vector<DressedLepton> muons     = apply<DressedLeptons>(event, "muons").dressedLeptons();
      const Jets& jets  = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      const Jets& fjets  = apply<FastJets>(event, "fjets").jetsByPt(Cuts::pT > 200*GeV && Cuts::abseta < 2.0);
      PseudoJets ljets;
      for (const Jet& fjet : fjets) { ljets += _trimmer(fjet); }
      sort(ljets.begin(), ljets.end(), [](PseudoJet const &l, PseudoJet const &r) { return l.pt() > r.pt(); });
      const FourMomentum& met = apply<MissingMomentum>(event, "MissingMomentum").missingMomentum();

      // Overlap removal
      for (const Jet& jet : jets) {
        ifilter_discard(electrons, deltaRLess(jet, 0.4, RAPIDITY));
        ifilter_discard(muons, deltaRLess(jet, 0.4, RAPIDITY));
      }

      // Reconstruct event

      // Lepton l
      size_t n_el_25 = 0, n_el_27 = 0;
      for (const DressedLepton& electron : electrons ) {
        if (electron.pT() >= 25*GeV)  ++n_el_25;
        if (electron.pT() >= 27*GeV)  ++n_el_27;
      }
      size_t n_mu_25 = 0, n_mu_27 = 0;
      for (const DressedLepton& muon : muons ) {
        if (muon.pT() >= 25*GeV)  ++n_mu_25;
        if (muon.pT() >= 27*GeV)  ++n_mu_27;
      }
      if ((n_el_25 + n_mu_25 != 1) || (n_el_27 + n_mu_27 != 1))  vetoEvent;
      DressedLepton lepton = (n_el_27 == 1)? electrons[0] : muons[0];
      FourMomentum l = lepton.mom();
      int lep_charge = lepton.charge();

      // Neutrino nu
      FourMomentum nu = getNeutrino(l, met);

      // Hadronic top candidate jh
      int jh_idx = -1;
      for (size_t ijet = 0; ijet < ljets.size(); ++ijet) {
        Jet ljet = Jet(ljets[ijet]);
        if (ljet.pT() < 350*GeV)     continue;
        if (ljet.abseta() > 2.0)     continue;
        if (ljet.mass() < 140*GeV)   continue;
        if (deltaPhi(ljet,l) < 1.0)  continue;
        bool btagged = false;
        for (const Jet& jet : jets) {
          if ( jet.bTagged(Cuts::pT > 5*GeV) ) {
            if ( deltaR(ljet, jet) < 1.0 ) {
              btagged = true;
              break;
            }
          }
        }
        if (btagged) {
          jh_idx = ijet;
          break;
        }
      }
      if ( jh_idx == -1 ) vetoEvent;
      FourMomentum jh = Jet(ljets[jh_idx]).mom();

      // Leptonic top b-jet candidate jl
      int jl_idx = -1;
      for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
        if ( !jets[ijet].bTagged(Cuts::pT > 5*GeV) ) continue;
        if ( deltaR(jets[ijet],  l) > 2.0 ) continue;
        if ( deltaR(jets[ijet], jh) < 1.5 ) continue;
        jl_idx = ijet;
        break;
      }
      if ( jl_idx == -1) {
        for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
          if ( deltaR(jets[ijet],  l) > 2.0 ) continue;
          if ( deltaR(jets[ijet], jh) < 1.5 ) continue;
          jl_idx = ijet;
          break;
        }
      }
      if  ( jl_idx == -1 ) vetoEvent;
      FourMomentum jl = Jet(jets[jl_idx]).mom();

      // b-tagging
      size_t n_btagged = 0;
      size_t n_btagged_matched = 1; // Large-jet jh is b-tagged
      for (const Jet& jet : jets ) {
        if (jet.bTagged(Cuts::pT > 5*GeV))  ++n_btagged;
      }
      if ( jets[jl_idx].bTagged(Cuts::pT > 5*GeV) )  ++n_btagged_matched;
      if ( n_btagged >= 2 && n_btagged_matched < 2 )  vetoEvent;

      // Associated jet candidate ja
      int ja_idx = -1;
      for (int ijet = 0; ijet < int(jets.size()); ++ijet) {
        if ( ijet == jl_idx ) continue;
        if ( jets[ijet].pT() < 100*GeV ) continue;
        if ( deltaR(jets[ijet], jh) < 1.5 ) continue;
        if ( deltaR(jets[ijet],  l) < 0.4 ) continue;
        ja_idx = ijet;
        break;
      }
      if ( ja_idx == -1 ) vetoEvent;
      FourMomentum ja = jets[ja_idx].mom();

      FourMomentum thad = jh;
      FourMomentum tlep = l+nu+jl;
      FourMomentum top = lep_charge > 0 ? tlep : thad;
      FourMomentum tbar = lep_charge > 0 ? thad : tlep;
      FourMomentum ttbar = top + tbar;
      FourMomentum ttj = top + tbar + ja;

      // Boost into ttj reference frame
      FourMomentum ttj_inv( ttj.E(), -ttj.px(), -ttj.py(), -ttj.pz() );
      Vector3 boostVector = ttj_inv.betaVec();
      LorentzTransform lt_boost;
      lt_boost.setBetaVec( boostVector );
      FourMomentum top_boosted = lt_boost.transform( top );
      FourMomentum tbar_boosted = lt_boost.transform( tbar );
      FourMomentum ja_boosted = lt_boost.transform( ja );

      // Get observables
      const double deltaE = top_boosted.E() - tbar_boosted.E();
      const double thetaj_opt = ttj.rapidity() > 0 ? ja_boosted.theta() : pi - ja_boosted.theta();

      // Fill auxiliary histograms
      _h[deltaE > 0 ? "pos" : "neg"]->fill(thetaj_opt);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h, crossSection()/picobarn/sumW());

      // Calculate differential energy asymmetry
      asymm(_h["pos"], _h["neg"], _asymm);

    }

    /// @}


  private:

    fastjet::Filter _trimmer;

    // Histograms
    map<string, Histo1DPtr> _h;
    Scatter2DPtr _asymm;


    static double delta2_fcn(const MendelMin::Params& p, const MendelMin::Params& pfix) {
      double delta2 = 0;
      double alpha = p[0]*6.30-3.15; // Map p[0] in [0,1] to alpha in [-3.15,3.15]
      double r = pfix[0];
      double dphi = pfix[1];
      double l_pt = pfix[2];
      double l_m = pfix[3];
      double n_px = pfix[4];
      double n_py = pfix[5];
      r /= sqrt(l_pt * l_pt + l_m * l_m) - l_pt * cos(dphi + alpha);
      FourMomentum neut(0.0, n_px, n_py, 0.0); // E, px, py, pz
      neut.setE(neut.p());
      FourMomentum neut_new( 0.0, r * neut.p() * cos(neut.phi() + alpha), r * neut.p() * sin(neut.phi() + alpha), 0.0 );
      neut_new.setE(neut_new.p());
      delta2 = pow((neut_new.px() - neut.px()), 2) + pow((neut_new.py() - neut.py()), 2);
      return delta2;
    }


    FourMomentum getNeutrino(const FourMomentum& lepton, const FourMomentum& met) {
      const double m_mWpdg = 80.4*GeV;
      double pxNu = met.px();
      double pyNu = met.py();
      double ptNu = met.pt();
      double pzNu;

      double c1 = pow(m_mWpdg,2) - pow(lepton.mass(),2) + 2 * (lepton.px() * pxNu + lepton.py() * pyNu);
      double b1 = 2 * lepton.pz();
      double A = 4*pow(lepton.E(),2) - b1*b1;
      double B = -2 * c1 * b1;
      double C = 4 * pow(lepton.E(), 2) * ptNu * ptNu - c1 * c1;
      double discr = B*B - 4*A*C;
      double r = 1;
      double sol1, sol2;
      if (discr > 0){
        sol1 = (-B + sqrt(discr)) / (2*A);
        sol2 = (-B - sqrt(discr)) / (2*A);
      }
      else {
        // fitAlpha
        std::valarray<double> pfix  = { (m_mWpdg * m_mWpdg - lepton.mass() * lepton.mass()) / (2 * ptNu), met.phi() - lepton.phi(),
                                        lepton.pt(), lepton.mass(), pxNu, pyNu };
        MendelMin mm(delta2_fcn, 1, pfix);
        mm.evolve(100);
        valarray<double> fittest = mm.fittest();

        const double alpha = fittest[0]*6.30-3.15; // map p[0] in [0,1] to alpha in [-3.15,3.15]
        const double dphi = met.phi() - lepton.phi();
        r  = ( pow(m_mWpdg,2) - pow(lepton.mass(),2) );
        r /= (2 * ptNu * (sqrt(pow(lepton.pt(),2) + pow(lepton.mass(),2)) - lepton.pt() * cos(dphi + alpha)));

        const double old_p = ptNu;
        const double old_phi = met.phi();
        pxNu = r * old_p * cos(old_phi + alpha);
        pyNu = r * old_p * sin(old_phi + alpha);
        ptNu = sqrt (pxNu*pxNu + pyNu*pyNu);

        c1 = pow(m_mWpdg,2) - pow(lepton.mass(),2) + 2 * (lepton.px() * pxNu + lepton.py() * pyNu);
        B = -2 * c1 * b1;
        C = 4 * pow(lepton.E(),2) * ptNu * ptNu - c1 * c1;
        discr = B*B - 4*A*C;

        sol1 = -B / (2*A);
        sol2 = -B / (2*A);
      }
      // useSmallestPz
      pzNu = (fabs(sol1) > fabs(sol2)) ? sol2 : sol1;

      FourMomentum nu( sqrt(sqr(pxNu) + sqr(pyNu) + sqr(pzNu)), pxNu, pyNu, pzNu);

      return nu;
    }

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2021_I1941095);

}
