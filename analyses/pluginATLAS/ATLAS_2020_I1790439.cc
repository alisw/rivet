// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"


namespace Rivet {

  /// @brief H->ZZ->4l 13 TeV analysis
  class ATLAS_2020_I1790439 : public Analysis {
  public:
    /// Default constructor
    ATLAS_2020_I1790439()
      : Analysis("ATLAS_2020_I1790439"),
        MGME()
    { }

    void init() {

      /// Dressed leptons
      Cut cut_lep = (Cuts::abseta < 2.7) && (Cuts::pT > 5*GeV);
      PromptFinalState prompt_photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState prompt_leptons(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      DressedLeptons dLeptons(prompt_photons, prompt_leptons, 0.1, cut_lep, true);
      declare(dLeptons, "AllLeptons");

      /// Jet inputs
      FinalState fs_jet(Cuts::abseta < 5.0);
      VetoedFinalState jet_input(fs_jet);

      // reject all leptons dressed with only prompt photons from jet input
      FinalState leptons(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      DressedLeptons reject_leptons(prompt_photons, leptons, 0.1, Cuts::open(), true);
      jet_input.addVetoOnThisFinalState(reject_leptons);

      // reject prompt invisibles, including from tau decays
      VetoedFinalState invis_fs_jet(fs_jet);
      invis_fs_jet.addVetoOnThisFinalState(VisibleFinalState(fs_jet));
      PromptFinalState invis_pfs_jet = PromptFinalState(invis_fs_jet, true);
      jet_input.addVetoOnThisFinalState(invis_pfs_jet);

      // declare jets
      FastJets jets(jet_input, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::DECAY);
      declare(jets, "Jets");

      // Book histograms
      book(_h["xs_flavor"],           3, 1, 1);
      book(_h["H4l_pt"],              5, 1, 1);
      book(_h["Z1_m"],                7, 1, 1);
      book(_h["Z2_m"],                9, 1, 1);
      book(_h["abshiggs_y"],          11, 1, 1);
      book(_h["abscthstr"],           13, 1, 1);
      book(_h["cth1"],                15, 1, 1);
      book(_h["cth2"],                17, 1, 1);
      book(_h["phi"],                 19, 1, 1);
      book(_h["phi1"],                21, 1, 1);
      book(_h["n_jets"],              23, 1, 1);
      book(_h["n_jets_incl"],         25, 1, 1);
      book(_h["n_bjets"],             26, 1, 1);
      book(_h["jet_pt_leading"],      28, 1, 1);
      book(_h["jet_pt_subleading"],   30, 1, 1);
      book(_h["dijet_m"],             32, 1, 1);
      book(_h["dijet_deltaeta"],      34, 1, 1);
      book(_h["dijet_deltaphi"],      36, 1, 1);
      book(_h["pt4lj"],               38, 1, 1);
      book(_h["pt4ljj"],              40, 1, 1);
      book(_h["m4lj"],                42, 1, 1);
      book(_h["m4ljj"],               44, 1, 1);
      book(_h["m12vsm34"],            46, 1, 1);
      book(_h["m12vsm34_2l2m"],       48, 1, 1);
      book(_h["m12vsm34_2l2e"],       49, 1, 1);
      book(_h["pt4lvy4l_0_0p5"],      51, 1, 1);
      book(_h["pt4lvy4l_0p5_1"],      51, 1, 3);
      book(_h["pt4lvy4l_1_1p5"],      51, 1, 5);
      book(_h["pt4lvy4l_1p5_2p5"],    51, 1, 7);
      book(_h["pt4lvnjet_0"],         53, 1, 1);
      book(_h["pt4lvnjet_1"],         53, 1, 3);
      book(_h["pt4lvnjet_2"],         53, 1, 5);
      book(_h["pt4lvnjet_3"],         53, 1, 7);
      book(_h["pt4lvpt4lj"],          55, 1, 1);
      book(_h["pt4ljvm4lj"],          57, 1, 1);
      book(_h["pt4lvptj0"],           59, 1, 1);
      book(_h["ptj0vyj0"],            61, 1, 1);
      book(_h["ptj0vptj1"],           63, 1, 1);
      book(_h["Z1_m_4l"],             65, 1, 1);
      book(_h["Z1_m_2l2l"],           66, 1, 1);
      book(_h["Z2_m_4l"],             68, 1, 1);
      book(_h["Z2_m_2l2l"],           69, 1, 1);
      book(_h["phi_4l"],              71, 1, 1);
      book(_h["phi_2l2l"],            72, 1, 1);
      book(_h["m12vsm34_4l"],         74, 1, 1);
      book(_h["m12vsm34_2l2l"],       75, 1, 1);
    }

    /// Do the per-event analysis
    void analyze(const Event& e) {

      _h["xs_flavor"]->fill(9);

      const std::vector<DressedLepton>& all_leps = apply<DressedLeptons>(e, "AllLeptons").dressedLeptons();
      unsigned int n_parts = all_leps.size();
      unsigned int n_OSSF_pairs = 0;
      std::vector<Zstate> dileptons;
      Particles leptons;

      // Form Z candidate (opposite-sign same-flavour) lepton pairs
      for (unsigned int i = 0; i < n_parts; i++) {
        for (unsigned int j = i + 1; j < n_parts; j++) {
          if (isOSSF(all_leps[i], all_leps[j])){
            n_OSSF_pairs += 1;
            // Set positive charge first for later ME calculation
            if (all_leps[i].charge() > 0) {
              dileptons.push_back(Zstate( ParticlePair(all_leps[i], all_leps[j]) ));
            } else {
              dileptons.push_back(Zstate( ParticlePair(all_leps[j], all_leps[i]) ));
            }
          }
        }
      }

      // At least two pairs required to select ZZ->llll final state
      if (n_OSSF_pairs < 2) vetoEvent;


      // Form the quadruplet of two lepon pairs passing kinematics cuts
      std::vector<Quadruplet> quadruplets;
      for (unsigned int i = 0; i < dileptons.size(); i++) {
        for (unsigned int j = i+1; j < dileptons.size(); j++) {
          // Only use unique leptons
          if (isSame( dileptons[i].first , dileptons[j].first  )) continue;
          if (isSame( dileptons[i].first , dileptons[j].second )) continue;
          if (isSame( dileptons[i].second, dileptons[j].first  )) continue;
          if (isSame( dileptons[i].second, dileptons[j].second )) continue;

          leptons.clear();
          leptons.push_back( dileptons[i].first  );
          leptons.push_back( dileptons[i].second );
          leptons.push_back( dileptons[j].first  );
          leptons.push_back( dileptons[j].second );
          leptons = sortByPt(leptons);

          // Apply kinematic cuts
          if ( leptons[0].pt() < 20*GeV) continue;
          if ( leptons[1].pt() < 15*GeV) continue;
          if ( leptons[2].pt() < 10*GeV) continue;

          // Form the quad with pair closest to Z pole first
          if (dileptons[i].Zdist() < dileptons[j].Zdist()) {
            quadruplets.push_back(Quadruplet(dileptons[i], dileptons[j]));
          } else {
            quadruplets.push_back(Quadruplet(dileptons[j], dileptons[i]));
          }
        }
      }

      // Veto if no quad passes kinematic selection
      unsigned int n_quads = quadruplets.size();
      if(n_quads == 0) vetoEvent;

      // To resolve ambiguities in lepton pairing order quads by channel priority first, then m12 - mz and m34 - mz
      // The first in every channel is considered nominal
      std::sort(quadruplets.begin(), quadruplets.end(),
                [](const Quadruplet & q1, const Quadruplet & q2) {
                  if (q1.ch_priority == q2.ch_priority) {
                    // if rarely, Z1 the same distance from the Z pole, compare Z2
                    if (fabs( q1.Z1().Zdist() - q2.Z1().Zdist() ) < 1.e-5)
                      return q1.Z2().Zdist() < q2.Z2().Zdist();
                    else
                      return q1.Z1().Zdist() < q2.Z1().Zdist();
                  } else
                    return q1.ch_priority < q2.ch_priority;
                });


      // Select the best quad
      Particles leptons_sel4l;
      Quadruplet quadSel;
      float MEHZZ_max  = -999.;
      int prevQuadType = -999.;
      bool atleastonequadpassed = false;
      bool isNominalQuad        = false;
      bool extraLep             = false;

      for(unsigned int iquad = 0; iquad < n_quads; iquad++) {

        // Veto event if nominal quad was not selected in 4 lepton case
        if (n_parts == 4 && iquad > 0) vetoEvent;

        int quadType = (int) quadruplets[iquad].type();

        if (quadType != prevQuadType) {
          isNominalQuad = true;
        } else {
          isNominalQuad = false;
        }
        prevQuadType = quadType;

        Quadruplet & quad = quadruplets[iquad];

        // Z invariant mass requirements
        if (!(inRange(quad.Z1().mom().mass(), 50*GeV, 106*GeV))) continue;
        if (!(inRange(quad.Z2().mom().mass(), 12*GeV, 115*GeV))) continue;

        // Lepton separation and J/Psi veto
        bool b_pass_leptonseparation = true;
        bool b_pass_jpsi = true;
        leptons_sel4l.clear();
        leptons_sel4l.push_back(quad.Z1().first);
        leptons_sel4l.push_back(quad.Z1().second);
        leptons_sel4l.push_back(quad.Z2().first);
        leptons_sel4l.push_back(quad.Z2().second);

        for (unsigned int i = 0; i < 4; i++) {
          for (unsigned int j = i+1; j < 4; j++) {
            if ( deltaR( leptons_sel4l[i], leptons_sel4l[j]) < 0.1) b_pass_leptonseparation = false;
            if ( isOSSF(leptons_sel4l[i], leptons_sel4l[j]) && (leptons_sel4l[i].mom() + leptons_sel4l[j].mom()).mass() <= 5.*GeV) b_pass_jpsi = false;
          }
        }
        if(b_pass_leptonseparation == false || b_pass_jpsi == false) continue;

        // Only consider the event if at least one nominal quadruplet passes cuts
        if (isNominalQuad) atleastonequadpassed = true;

        // Direct selection for case with only 4 prompt leptons
        if (n_parts == 4 ) {
          quadSel = quad;
          break;
        }

        // In cases with extra leptons meeting further requirements use max ME to select quad
        float MEHZZ;
        if (quadType == 0 || quadType == 1) MEHZZ = MGME.Compute(quadruplets[iquad]) / 2;
        else MEHZZ = MGME.Compute(quadruplets[iquad]);

        // Check leptons other than the ones from the channel's nominal quadruplet
        // and don't need to recheck extraLep for a second channel
        if(isNominalQuad && !extraLep) {
          for (const Particle &lep: all_leps) {
            bool lep_in_quad = false;
            bool lep_too_close = false;

            // pT requirement
            if (lep.pt() < 12*GeV) continue;

            // skip if lepton included in quad or not isolated from quad leptons
            for (unsigned int i = 0; i < 4; i++) {
              if (isSame(lep, leptons_sel4l[i])) lep_in_quad = true ;
              if (deltaR(lep, leptons_sel4l[i]) < 0.1) lep_too_close = true ;
            }
            if (lep_in_quad || lep_too_close) continue;

            extraLep = true;
            break;
          }

          if(!extraLep) {
            // In case of no suitable extra leptons select the quad directly
            quadSel = quad;
            break;
          }
        }

        // Use ME to select the quad when there are extra leptons
        if (MEHZZ > MEHZZ_max) {
          quadSel = quad;
          MEHZZ_max = MEHZZ;
        }
      }

      if(!atleastonequadpassed) vetoEvent;

      // Veto if quad not in Higgs mass window
      FourMomentum Higgs = quadSel.mom();
      double H4l_mass     = Higgs.mass();
      if (!(inRange(H4l_mass, 105.*GeV, 160.*GeV))) vetoEvent;

      // Higgs observables
      double H4l_pt       = Higgs.pt();
      double H4l_rapidity = Higgs.absrapidity();
      LorentzTransform HRF_boost = LorentzTransform::mkFrameTransformFromBeta(Higgs.betaVec());

      FourMomentum Z1_in_HRF = HRF_boost.transform( quadSel.Z1().mom() );
      FourMomentum Z2_in_HRF = HRF_boost.transform( quadSel.Z2().mom() );
      double H4l_costheta = fabs(cos( Z1_in_HRF.theta()));
      double H4l_m12      = quadSel.Z1().mom().mass();
      double H4l_m34      = quadSel.Z2().mom().mass();

      FourMomentum v11_HRF = HRF_boost.transform( quadSel.Z1().second.mom() );
      FourMomentum v12_HRF = HRF_boost.transform( quadSel.Z1().first.mom() );
      FourMomentum v21_HRF = HRF_boost.transform( quadSel.Z2().second.mom() );
      FourMomentum v22_HRF = HRF_boost.transform( quadSel.Z2().first.mom() );

      Vector3 v11 = v11_HRF.p3();
      Vector3 v12 = v12_HRF.p3();
      Vector3 v21 = v21_HRF.p3();
      Vector3 v22 = v22_HRF.p3();
      Vector3 nz(0, 0, 1);
      Vector3 qz1 = Z1_in_HRF.p3();
      Vector3 n1p = v11.cross(v12).unit();
      Vector3 n2p = v21.cross(v22).unit();
      Vector3 nscp = nz.cross(qz1).unit();

      double H4l_Phi  = qz1.dot(n1p.cross(n2p)) / fabs(qz1.dot(n1p.cross(n2p))) * acos(- n1p.dot(n2p));
      double H4l_Phi1 = qz1.dot(n1p.cross(nscp)) / fabs(qz1.dot(n1p.cross(nscp))) * acos( n1p.dot(nscp));

      LorentzTransform Z1RF_boost = LorentzTransform::mkFrameTransformFromBeta(Z1_in_HRF.betaVec());
      LorentzTransform Z2RF_boost = LorentzTransform::mkFrameTransformFromBeta(Z2_in_HRF.betaVec());

      FourMomentum Z1_in_Z2RF = Z2RF_boost.transform( Z1_in_HRF );
      FourMomentum Z2_in_Z1RF = Z1RF_boost.transform( Z2_in_HRF );

      Vector3 Z1_p3 = Z1_in_Z2RF.p3();
      Vector3 Z2_p3 = Z2_in_Z1RF.p3();

      FourMomentum n_Z1 = Z1RF_boost.transform( v11_HRF );
      FourMomentum n_Z2 = Z2RF_boost.transform( v21_HRF );

      // angle is negative with its Z
      double H4l_cth1 = - cos( n_Z1.p3().angle(Z2_p3));
      double H4l_cth2 = - cos( n_Z2.p3().angle(Z1_p3));


      // Jet observables
      Jets jets = apply<FastJets>(e, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);

      // discard jets which overlap leptons
      for (const Particle& l: all_leps) ifilter_discard(jets, deltaRLess(l, 0.1));

      // collect b-tagged jets
      const Jets b_jets_btagged = filter_select(jets, hasBTag(Cuts::pT>5*GeV));

      unsigned int n_bjets = b_jets_btagged.size();
      unsigned int n_jets = jets.size();
      double leading_jet_pt = (n_jets > 0 ? jets[0].pT() : 0);
      double leading_jet_y = (n_jets > 0 ? jets[0].absrapidity() : 0);
      double subleading_jet_pt = (n_jets > 1 ? jets[1].pT() : 0);

      FourMomentum Dijets = {0,0,0,0};
      if(n_jets > 1) Dijets = jets[0].mom() + jets[1].mom();
      double mjj = (n_jets > 1 ? Dijets.mass() : -999);

      double dphijj = -1;
      if(n_jets > 1){
        if(jets[0].eta() >  jets[1].eta()) dphijj = jets[0].phi() - jets[1].phi();
        else dphijj = jets[1].phi() - jets[0].phi();
        if(dphijj < 0) dphijj = dphijj + TWOPI;
      }

      double detajj = (n_jets > 1 ? fabs(jets[0].eta() -  jets[1].eta()): -1);

      FourMomentum m4lj, m4ljj;
      if (n_jets > 0) m4lj  = Higgs + jets[0];
      if (n_jets > 1) m4ljj = Higgs + Dijets;

      double H4l_m4lj = m4lj.mass();
      double H4l_m4ljj = m4ljj.mass();
      double H4l_pt4lj = m4lj.pt();
      double H4l_pt4ljj = m4ljj.pt();

      // Fill histograms
      // Branching ratios for 4l and 2l2l channels
      double BR_SF = 0.00013;
      double BR_OF = 0.000118;
      if (inRange(H4l_mass, 115.*GeV, 130.*GeV)){
        if (quadSel.type() == Quadruplet::FlavCombi::mm
            || quadSel.type() == Quadruplet::FlavCombi::ee ) {
          _h["xs_flavor"]->fill((int)quadSel.type()+1, BR_SF);
          _h["xs_flavor"]->fill(5, BR_SF);
        }
        else if (quadSel.type() == Quadruplet::FlavCombi::em
                 || quadSel.type() == Quadruplet::FlavCombi::me ) {
          _h["xs_flavor"]->fill((int)quadSel.type()+1, BR_OF);
          _h["xs_flavor"]->fill(6, BR_OF);
        }
        _h["xs_flavor"]->fill(7, Br);
        _h["xs_flavor"]->fill(8, Br);
      }

      // Higgs variables
      for(const auto & p: std::map<std::string, double>{
          {"H4l_pt",     H4l_pt},
          {"Z1_m",       H4l_m12},
          {"Z2_m",       H4l_m34},
          {"abshiggs_y", H4l_rapidity},
          {"abscthstr",  H4l_costheta},
          {"cth1",       H4l_cth1},
          {"cth2",       H4l_cth2},
          {"phi",        H4l_Phi},
          {"phi1",       H4l_Phi1}})
        {
          _h[ p.first ]->fill(p.second);
        }

      // Jet variables
      if (n_jets <= 2) _h[ "n_jets" ]->fill(n_jets);
      else _h[ "n_jets" ]->fill(3.0);

      _h[ "n_jets_incl" ]->fill(0.0);
      if (n_jets >= 1) _h[ "n_jets_incl" ]->fill(1.0);
      if (n_jets >= 2) _h[ "n_jets_incl" ]->fill(2.0);
      if (n_jets >= 3) _h[ "n_jets_incl" ]->fill(3.0);

      if (n_jets == 0) _h[ "n_bjets" ]->fill(1.0);
      else if (n_bjets == 0) _h[ "n_bjets" ]->fill(2.0);
      else if (n_bjets >= 1) _h[ "n_bjets" ]->fill(3.0);

      if (n_jets == 0) {
        _h[ "jet_pt_leading" ]->fill(29.5);
        _h[ "pt4lj" ]->fill(-0.5);
        _h[ "m4lj" ]->fill(119.5);
      }
      else if(n_jets >= 1){
        _h[ "jet_pt_leading" ]->fill(leading_jet_pt);
        _h[ "pt4lj" ]->fill(H4l_pt4lj);
        _h[ "m4lj" ]->fill(H4l_m4lj);
      }

      if (n_jets < 2) {
        _h[ "jet_pt_subleading" ]->fill(29.5);
        _h[ "dijet_m" ]->fill(-0.5);
        _h[ "dijet_deltaeta" ]->fill(-0.5);
        _h[ "dijet_deltaphi" ]->fill(-0.5);
        _h[ "pt4ljj" ]->fill(-0.5);
        _h[ "m4ljj" ]->fill(179.5);
      }
      else if (n_jets >= 2) {
        _h[ "jet_pt_subleading" ]->fill(subleading_jet_pt);
        _h[ "dijet_m" ]->fill(mjj);
        _h[ "dijet_deltaeta" ]->fill(detajj);
        _h[ "dijet_deltaphi" ]->fill(dphijj);
        _h[ "pt4ljj" ]->fill(H4l_pt4ljj);
        _h[ "m4ljj" ]->fill(H4l_m4ljj);
      }

        // m12 vs m34 (all channels)
        if(H4l_m12 < 82 && H4l_m34 < 32) _h["m12vsm34"]->fill(1.);
        else if(H4l_m12 < 74 && H4l_m34 > 32) _h["m12vsm34"]->fill(2.);
        else if(H4l_m12 > 74 && H4l_m34 > 32) _h["m12vsm34"]->fill(3.);
        else if(H4l_m12 > 82 && H4l_m34 < 32 && H4l_m34 > 24) _h["m12vsm34"]->fill(4.);
        else if(H4l_m12 > 82 && H4l_m34 < 24) _h["m12vsm34"]->fill(5.);

        if (quadSel.type() == Quadruplet::FlavCombi::em
            || quadSel.type() == Quadruplet::FlavCombi::mm ){
          if(H4l_m12 < 82 && H4l_m34 < 32) _h["m12vsm34_2l2m"]->fill(1.);
          else if(H4l_m12 < 74 && H4l_m34 > 32) _h["m12vsm34_2l2m"]->fill(2.);
          else if(H4l_m12 > 74 && H4l_m34 > 32) _h["m12vsm34_2l2m"]->fill(3.);
          else if(H4l_m12 > 82 && H4l_m34 < 32 && H4l_m34 > 24) _h["m12vsm34_2l2m"]->fill(4.);
          else if(H4l_m12 > 82 && H4l_m34 < 24) _h["m12vsm34_2l2m"]->fill(5.);
        }
        else if (quadSel.type() == Quadruplet::FlavCombi::me
                 || quadSel.type() == Quadruplet::FlavCombi::ee ){
          if(H4l_m12 < 82 && H4l_m34 < 32) _h["m12vsm34_2l2e"]->fill(1.);
          else if(H4l_m12 < 74 && H4l_m34 > 32) _h["m12vsm34_2l2e"]->fill(2.);
          else if(H4l_m12 > 74 && H4l_m34 > 32) _h["m12vsm34_2l2e"]->fill(3.);
          else if(H4l_m12 > 82 && H4l_m34 < 32 && H4l_m34 > 24) _h["m12vsm34_2l2e"]->fill(4.);
          else if(H4l_m12 > 82 && H4l_m34 < 24) _h["m12vsm34_2l2e"]->fill(5.);
        }

        // m12 vs m34 (4l channels only)
        if (quadSel.type() == Quadruplet::FlavCombi::ee
            || quadSel.type() == Quadruplet::FlavCombi::mm ){
          if(H4l_m12 < 82 && H4l_m34 < 32) _h["m12vsm34_4l"]->fill(1.);
          else if(H4l_m12 < 74 && H4l_m34 > 32) _h["m12vsm34_4l"]->fill(2.);
          else if(H4l_m12 > 74 && H4l_m34 > 32) _h["m12vsm34_4l"]->fill(3.);
          else if(H4l_m12 > 82 && H4l_m34 < 32 && H4l_m34 > 24) _h["m12vsm34_4l"]->fill(4.);
          else if(H4l_m12 > 82 && H4l_m34 < 24) _h["m12vsm34_4l"]->fill(5.);
          _h["Z1_m_4l"]->fill(H4l_m12);
          _h["Z2_m_4l"]->fill(H4l_m34);
          _h["phi_4l"]->fill(H4l_Phi);
        }

        // m12 vs m34 (2l2l channels only)
        if (quadSel.type() == Quadruplet::FlavCombi::me
            || quadSel.type() == Quadruplet::FlavCombi::em ){
          if(H4l_m12 < 82 && H4l_m34 < 32) _h["m12vsm34_2l2l"]->fill(1.);
          else if(H4l_m12 < 74 && H4l_m34 > 32) _h["m12vsm34_2l2l"]->fill(2.);
          else if(H4l_m12 > 74 && H4l_m34 > 32) _h["m12vsm34_2l2l"]->fill(3.);
          else if(H4l_m12 > 82 && H4l_m34 < 32 && H4l_m34 > 24) _h["m12vsm34_2l2l"]->fill(4.);
          else if(H4l_m12 > 82 && H4l_m34 < 24) _h["m12vsm34_2l2l"]->fill(5.);
          _h["Z1_m_2l2l"]->fill(H4l_m12);
          _h["Z2_m_2l2l"]->fill(H4l_m34);
          _h["phi_2l2l"]->fill(H4l_Phi);
        }

        // 2d differential variables
        if (0 < H4l_rapidity && H4l_rapidity < 0.5) _h["pt4lvy4l_0_0p5"]->fill(H4l_pt);
        else if (0.5 < H4l_rapidity && H4l_rapidity < 1) _h["pt4lvy4l_0p5_1"]->fill(H4l_pt);
        else if (1 < H4l_rapidity && H4l_rapidity < 1.5) _h["pt4lvy4l_1_1p5"]->fill(H4l_pt);
        else if (1.5 < H4l_rapidity && H4l_rapidity < 2.5) _h["pt4lvy4l_1p5_2p5"]->fill(H4l_pt);

        if (n_jets == 0) _h["pt4lvnjet_0"]->fill(H4l_pt);
        else if (n_jets == 1) _h["pt4lvnjet_1"]->fill(H4l_pt);
        else if (n_jets == 2) _h["pt4lvnjet_2"]->fill(H4l_pt);
        else if (n_jets > 2) _h["pt4lvnjet_3"]->fill(H4l_pt);

        if (n_jets == 0) {
          _h["pt4lvptj0"]->fill(1.);
          _h["pt4lvpt4lj"]->fill(1.);
          _h["pt4ljvm4lj"]->fill(1.);
          _h["ptj0vptj1"]->fill(1.);
          _h["ptj0vyj0"]->fill(1.);
        } else {

          if (     0  < H4l_pt4lj && H4l_pt4lj < 60  && 0   < H4l_pt && H4l_pt < 120)  _h["pt4lvpt4lj"]->fill(2.);
          else if (0  < H4l_pt4lj && H4l_pt4lj < 60  && 120 < H4l_pt && H4l_pt < 350)  _h["pt4lvpt4lj"]->fill(3.);
          else if (60 < H4l_pt4lj && H4l_pt4lj < 350 && 0   < H4l_pt && H4l_pt < 120)  _h["pt4lvpt4lj"]->fill(4.);
          else if (60 < H4l_pt4lj && H4l_pt4lj < 350 && 120 < H4l_pt && H4l_pt < 350)  _h["pt4lvpt4lj"]->fill(5.);

          if (     120 < H4l_m4lj && H4l_m4lj < 220  && 0   < H4l_pt4lj && H4l_pt4lj < 350) _h["pt4ljvm4lj"]->fill(2.);
          else if (220 < H4l_m4lj && H4l_m4lj < 350  && 0   < H4l_pt4lj && H4l_pt4lj < 60)  _h["pt4ljvm4lj"]->fill(3.);
          else if (220 < H4l_m4lj && H4l_m4lj < 350  && 60  < H4l_pt4lj && H4l_pt4lj < 350) _h["pt4ljvm4lj"]->fill(4.);
          else if (350 < H4l_m4lj && H4l_m4lj < 2000 && 0   < H4l_pt4lj && H4l_pt4lj < 350) _h["pt4ljvm4lj"]->fill(5.);

          if (     30 < leading_jet_pt  && leading_jet_pt < 60  && 0   < H4l_pt && H4l_pt < 80)  _h["pt4lvptj0"]->fill(2.);
          else if (30 < leading_jet_pt  && leading_jet_pt < 60  && 80  < H4l_pt && H4l_pt < 350) _h["pt4lvptj0"]->fill(3.);
          else if (60 < leading_jet_pt  && leading_jet_pt < 120 && 0   < H4l_pt && H4l_pt < 120) _h["pt4lvptj0"]->fill(4.);
          else if (60 < leading_jet_pt  && leading_jet_pt < 120 && 120 < H4l_pt && H4l_pt < 350) _h["pt4lvptj0"]->fill(5.);
          else if (120 < leading_jet_pt && leading_jet_pt < 350 && 0   < H4l_pt && H4l_pt < 120) _h["pt4lvptj0"]->fill(6.);
          else if (120 < leading_jet_pt && leading_jet_pt < 350 && 120 < H4l_pt && H4l_pt < 350) _h["pt4lvptj0"]->fill(7.);

          if (     30 < leading_jet_pt && leading_jet_pt < 120 &&  0   < leading_jet_y && leading_jet_y < 0.8) _h["ptj0vyj0"]->fill(2.);
          else if (30 < leading_jet_pt && leading_jet_pt < 120 &&  0.8 < leading_jet_y && leading_jet_y < 1.7) _h["ptj0vyj0"]->fill(3.);
          else if (30 < leading_jet_pt && leading_jet_pt < 120 &&  1.7 < leading_jet_y                       ) _h["ptj0vyj0"]->fill(4.);
          else if (120 < leading_jet_pt && leading_jet_pt < 350 && 0   < leading_jet_y && leading_jet_y < 1.7) _h["ptj0vyj0"]->fill(5.);
          else if (120 < leading_jet_pt && leading_jet_pt < 350 && 1.7 < leading_jet_y                       ) _h["ptj0vyj0"]->fill(6.);

          if (     n_jets == 1 && 30 < leading_jet_pt && leading_jet_pt < 60)  _h["ptj0vptj1"]->fill(2.);
          else if (n_jets == 1 && 60 < leading_jet_pt && leading_jet_pt < 350) _h["ptj0vptj1"]->fill(3.);
          else if (30 < leading_jet_pt && leading_jet_pt < 60  && 30 < subleading_jet_pt && subleading_jet_pt < 60)  _h["ptj0vptj1"]->fill(4.);
          else if (60 < leading_jet_pt && leading_jet_pt < 350 && 30 < subleading_jet_pt && subleading_jet_pt < 60)  _h["ptj0vptj1"]->fill(5.);
          else if (60 < leading_jet_pt && leading_jet_pt < 350 && 60 < subleading_jet_pt && subleading_jet_pt < 350) _h["ptj0vptj1"]->fill(6.);
        }
    }


    void finalize() {

      const double sf = crossSection() / femtobarn  / sumOfWeights();
      for (auto hist : _h) {
        if( hist.first == "xs_flavor"){
          scale(hist.second, sf);
        } else {
          scale(hist.second, sf * Br);
        }

        // Scale individual bins which have been widened for visability
        if (hist.first == "jet_pt_leading" || hist.first == "jet_pt_subleading") {
          hist.second->bin(0).scaleW(30);
        }
        else if (hist.first == "dijet_m") {
          hist.second->bin(0).scaleW(500);
        }
        else if (hist.first == "pt4lj" || hist.first == "pt4ljj") {
          hist.second->bin(0).scaleW(60);
        }
        else if (hist.first == "m4lj") {
          hist.second->bin(0).scaleW(120);
        }
        else if (hist.first == "m4ljj") {
          hist.second->bin(0).scaleW(180);
        }

      }
    }

  private:

    // Br(H-->ZZ) * BR(ZZ-->4l)
    const double Br = 0.02641 *  0.004736842;

    map<string, Histo1DPtr> _h;

    /// Generic Z candidate
    struct Zstate : public ParticlePair {
      Zstate() { }
      Zstate(ParticlePair _particlepair) : ParticlePair(_particlepair) { }
      FourMomentum mom() const { return first.momentum() + second.momentum(); }
      double Zdist() const { return fabs(mom().mass() -  91.1876*GeV); }
      int flavour() const { return first.abspid(); }
    };

    /// Generic quadruplet
    struct Quadruplet {

      // find out which type it is: 4mu = 0, 4e = 1, 2mu2e = 2, 2e2mu = 3 (mm, ee, me, em)
      // channel priority is 4m, 2e2m, 2m2e, 4e
      enum class FlavCombi { mm=0, ee, me, em, undefined };

      Quadruplet() { }
      Quadruplet(Zstate z1, Zstate z2) : _z1(z1), _z2(z2) {
        if (     _z1.flavour() == 13 && _z2.flavour() == 13) { _type = FlavCombi::mm; ch_priority = 0;}
        else if (_z1.flavour() == 11 && _z2.flavour() == 11) { _type = FlavCombi::ee; ch_priority = 3;}
        else if (_z1.flavour() == 13 && _z2.flavour() == 11) { _type = FlavCombi::me; ch_priority = 2;}
        else if (_z1.flavour() == 11 && _z2.flavour() == 13) { _type = FlavCombi::em; ch_priority = 1;}
        else  {_type = FlavCombi::undefined;}
      }

      Quadruplet(Quadruplet const & quad) :
        _z1(quad._z1),
        _z2(quad._z2),
        _type(quad._type),
        ch_priority(quad.ch_priority) {}

      Zstate _z1, _z2;
      FlavCombi _type;
      int ch_priority;

      const Zstate& Z1() const { return _z1; }
      const Zstate& Z2() const { return _z2; }
      FourMomentum mom() const { return _z1.mom() + _z2.mom(); }
      FlavCombi type() const {return _type; }
    };


    // save and calculate parameters
    class Parameters_heft {

    public:

      // Model parameters independent of aS
      double mdl_WH, mdl_WZ,
        aS, mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB,
        mdl_MP, mdl_conjg__CKM3x3, mdl_CKM3x3, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_MH__exp__4, mdl_MT__exp__4, mdl_MH__exp__2,
        mdl_MT__exp__2,
        mdl_MH__exp__6, mdl_MT__exp__6, mdl_aEW, mdl_MW, mdl_ee,
        mdl_MW__exp__2, mdl_sw2, mdl_cw,
        mdl_sw,
        mdl_v, mdl_ee__exp__2;
      std::complex<double> mdl_complexi;
      // Model parameters dependent on aS
      double mdl_GH ;
      // Model couplings independent of aS
      std::complex<double> GC_40, GC_54, GC_73;
      // Model couplings dependent on aS
      std::complex<double> GC_13;

      // Set parameters and couplings that are unchanged during the run
      Parameters_heft() {
        mdl_WH =  6.382339e-03;
        mdl_WZ =  2.441404e+00;
        aS =  1.180000e-01;
        mdl_Gf =  1.166390e-05;
        aEWM1 =  1.325070e+02;
        mdl_MH = 1.250000e+02;
        mdl_MZ = 9.118800e+01;
        mdl_MT = 1.730000e+02;
        mdl_complexi = std::complex<double> (0., 1.);
        mdl_MZ__exp__2 = pow(mdl_MZ, 2.);
        mdl_MZ__exp__4 = pow(mdl_MZ, 4.);
        mdl_MH__exp__2 = pow(mdl_MH, 2.);
        mdl_MT__exp__4 = pow(mdl_MT, 4.);
        mdl_MH__exp__4 = pow(mdl_MH, 4.);
        mdl_MT__exp__2 = pow(mdl_MT, 2.);
        mdl_MH__exp__6 = pow(mdl_MH, 6.);
        mdl_MT__exp__6 = pow(mdl_MT, 6.);
        mdl_aEW = 1./aEWM1;
        mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI * mdl_MZ__exp__2)/(mdl_Gf * sqrt(2.))));                                     mdl_ee = 2. * sqrt(mdl_aEW) * sqrt(M_PI);
        mdl_MW__exp__2 = pow(mdl_MW, 2.);
        mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2;
        mdl_cw = sqrt(1. - mdl_sw2);
        mdl_sw = sqrt(mdl_sw2);
        mdl_v = (2. * mdl_MW * mdl_sw)/mdl_ee;
        mdl_ee__exp__2 = pow(mdl_ee, 2.);
        GC_40 = -(mdl_ee * mdl_complexi * mdl_cw)/(2. * mdl_sw);
        GC_54 = (mdl_ee * mdl_complexi * mdl_sw)/(2. * mdl_cw);
        GC_73 = mdl_ee__exp__2 * mdl_complexi * mdl_v + ((1. - mdl_sw2) * mdl_ee__exp__2 * mdl_complexi * mdl_v)/(2. * mdl_sw2) +
          (mdl_ee__exp__2 * mdl_complexi * mdl_sw2 * mdl_v)/(2. * (1. - mdl_sw2));
      }

      // Set Mass
      void set4lepMass(double m_m4l){
        mdl_MH = m_m4l;
        mdl_WH = setWidth(m_m4l);
        mdl_MH__exp__2 = pow(mdl_MH, 2.);
        mdl_MH__exp__4 = pow(mdl_MH, 4.);
        mdl_MH__exp__6 = pow(mdl_MH, 6.);
        mdl_GH = -(4 * aS * M_PI * (1. + (13. * mdl_MH__exp__6)/(16800. * mdl_MT__exp__6) + mdl_MH__exp__4/(168. * mdl_MT__exp__4) +
                                    (7. * mdl_MH__exp__2)/(120. * mdl_MT__exp__2)))/(12. * pow(M_PI, 2.) * mdl_v);
        GC_13 = -(mdl_complexi * mdl_GH);
      }

    private:

      // Set Width
      long double setWidth(double m_m4l){
        long double Higgs_width_Poly_Fit_Zone1_coeff0  = -1.450308902710193E+03;
        long double Higgs_width_Poly_Fit_Zone1_coeff1  =  1.129291251156317E+02;
        long double Higgs_width_Poly_Fit_Zone1_coeff2  = -3.893063071316150E+00;
        long double Higgs_width_Poly_Fit_Zone1_coeff3  =  7.798666884832531E-02;
        long double Higgs_width_Poly_Fit_Zone1_coeff4  = -1.000455877406390E-03;
        long double Higgs_width_Poly_Fit_Zone1_coeff5  =  8.523735379647125E-06;
        long double Higgs_width_Poly_Fit_Zone1_coeff6  = -4.823164754652171E-08;
        long double Higgs_width_Poly_Fit_Zone1_coeff7  =  1.747954506786346E-10;
        long double Higgs_width_Poly_Fit_Zone1_coeff8  = -3.681723572169337E-13;
        long double Higgs_width_Poly_Fit_Zone1_coeff9  =  3.434207075968898E-16;

        long double Higgs_width_Poly_Fit_Zone2_coeff0  =  2.563291882845993E+02;
        long double Higgs_width_Poly_Fit_Zone2_coeff1  = -1.037082025855304E+01;
        long double Higgs_width_Poly_Fit_Zone2_coeff2  =  1.780260502696301E-01;
        long double Higgs_width_Poly_Fit_Zone2_coeff3  = -1.720311784419889E-03;
        long double Higgs_width_Poly_Fit_Zone2_coeff4  =  1.038418605369741E-05;
        long double Higgs_width_Poly_Fit_Zone2_coeff5  = -4.092496883922424E-08;
        long double Higgs_width_Poly_Fit_Zone2_coeff6  =  1.067667966800388E-10;
        long double Higgs_width_Poly_Fit_Zone2_coeff7  = -1.823343280081685E-13;
        long double Higgs_width_Poly_Fit_Zone2_coeff8  =  1.955637395597351E-16;
        long double Higgs_width_Poly_Fit_Zone2_coeff9  = -1.193287048560413E-19;
        long double Higgs_width_Poly_Fit_Zone2_coeff10 =  3.156196649452213E-23;

        long double Higgs_width_Poly_Fit_Zone3_coeff0  = -5.255605465437446E+02;
        long double Higgs_width_Poly_Fit_Zone3_coeff1  =  1.036972988796150E+01;
        long double Higgs_width_Poly_Fit_Zone3_coeff2  = -6.817022987365029E-02;
        long double Higgs_width_Poly_Fit_Zone3_coeff3  =  1.493275723660056E-04;

        long double m_m4l__2 = m_m4l * m_m4l;
        long double m_m4l__3 = m_m4l__2 * m_m4l;

        if( m_m4l < 156.5 ) return ( Higgs_width_Poly_Fit_Zone1_coeff0
                + Higgs_width_Poly_Fit_Zone1_coeff1 * m_m4l
                + Higgs_width_Poly_Fit_Zone1_coeff2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone1_coeff3 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone1_coeff4 * m_m4l__2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone1_coeff5 * m_m4l__2 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone1_coeff6 * m_m4l__2 * m_m4l__2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone1_coeff7 * m_m4l__2 * m_m4l__2 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone1_coeff8 * m_m4l__2 * m_m4l__2 * m_m4l__2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone1_coeff9 * m_m4l__2 * m_m4l__2 * m_m4l__2 * m_m4l__3 );
            else if( m_m4l >= 156.5 && m_m4l <= 162 ) return ( Higgs_width_Poly_Fit_Zone3_coeff0
                + Higgs_width_Poly_Fit_Zone3_coeff1 * m_m4l
                + Higgs_width_Poly_Fit_Zone3_coeff2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone3_coeff3 * m_m4l__3 );
            else return ( Higgs_width_Poly_Fit_Zone2_coeff0
                + Higgs_width_Poly_Fit_Zone2_coeff1 * m_m4l
                + Higgs_width_Poly_Fit_Zone2_coeff2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone2_coeff3 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone2_coeff4 * m_m4l__2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone2_coeff5 * m_m4l__2 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone2_coeff6 * m_m4l__2 * m_m4l__2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone2_coeff7 * m_m4l__2 * m_m4l__2 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone2_coeff8 * m_m4l__2 * m_m4l__2 * m_m4l__2 * m_m4l__2
                + Higgs_width_Poly_Fit_Zone2_coeff9 * m_m4l__2 * m_m4l__2 * m_m4l__2 * m_m4l__3
                + Higgs_width_Poly_Fit_Zone2_coeff10 * m_m4l__2 * m_m4l__2 * m_m4l__2 * m_m4l__2 * m_m4l__2 );
      }

    };  // class Parameters_heft

    // calculate LO ME
    class CPPProcess_P0_Sigma_heft_pp_H_ZZ_4l_heft_gg_epemmupmum {

    public:

      // Constructor.
      CPPProcess_P0_Sigma_heft_pp_H_ZZ_4l_heft_gg_epemmupmum() :
        pars(),
        pout(4, vector<double>(4, 0.)) { }

      // Destructor.
      virtual ~CPPProcess_P0_Sigma_heft_pp_H_ZZ_4l_heft_gg_epemmupmum() {}

      float Compute(const Quadruplet& quad) {
        FourMomentum cms = quad.mom();
        // use the cms mass as a value for MH
        pars.set4lepMass(cms.mass());

        const FourMomentum* fermionsMom[4] = {
          &quad.Z1().first.momentum(),
          &quad.Z1().second.momentum(),
          &quad.Z2().first.momentum(),
          &quad.Z2().second.momentum()
        };

        // boost to center-of-mass frame
        LorentzTransform HRF_boost = LorentzTransform::mkFrameTransformFromBeta(cms.betaVec());
        for(std::size_t i = 0; i < 4; ++i) {
          FourMomentum tmpMom = HRF_boost.transform(*fermionsMom[i]);
          pout[i][0] = tmpMom.E();
          pout[i][1] = tmpMom.px();
          pout[i][2] = tmpMom.py();
          pout[i][3] = tmpMom.pz();
        }

        // Evaluate matrix element
        return sigmaKin();
      }

    private:

      // Calculate flavour-independent parts of cross section. Evaluate |M|^2, part independent of
      // incoming flavour. Return matrix element

      const double sigmaKin() {
        // Local variables and constants
        static const int ncomb = 16;
        // Helicities for the process
        static const int helicities[ncomb][4] = {
          {-1, -1, -1, -1}, {-1, -1, -1,  1}, {-1, -1,  1, -1}, {-1, -1,  1,  1},
          {-1,  1, -1, -1}, {-1,  1, -1,  1}, {-1,  1,  1, -1}, {-1,  1,  1,  1},
          { 1, -1, -1, -1}, { 1, -1, -1,  1}, { 1, -1,  1, -1}, { 1, -1,  1,  1},
          { 1,  1, -1, -1}, { 1,  1, -1,  1}, { 1,  1,  1, -1}, { 1,  1,  1,  1}
        };

        // Reset the matrix elements
        double matrix_element = 0.;

        // Calculate the matrix element for all helicities
        for(int ihel = 0; ihel < ncomb; ihel++) {
          matrix_element += matrix_gg_h_h_zz_z_epem_z_mupmum(helicities[ihel]);
        }

        // Denominators: spins, colors and identical particles
        return matrix_element /= 128.;
      }

      // Calculate wavefunctions and matrix elements for all subprocesses
      double matrix_gg_h_h_zz_z_epem_z_mupmum(const int hel[]){
        // Calculate all wavefunctions
        ixx(pout[0], hel[0], w[2], true);
        ixx(pout[1], hel[1], w[3], false);
        FFV2_4_3(w[2], w[3], pars.GC_40, pars.GC_54,
                 pars.mdl_MZ, pars.mdl_WZ, w[4]);
        ixx(pout[2], hel[2], w[5], true);
        ixx(pout[3], hel[3], w[6], false);
        FFV2_4_3(w[5], w[6], pars.GC_40, pars.GC_54,
                 pars.mdl_MZ, pars.mdl_WZ, w[7]);
        VVS2_3(w[4], w[7], pars.GC_73, pars.mdl_MH, pars.mdl_WH, w[8]);

        // Calculate all amplitudes
        // Amplitude(s) for diagram number 0
        std::complex<double> amp = VVS3_0(pars.mdl_MH, w[8], pars.GC_13);

        // Calculate color flows
        // Sum and square the color flows to get the matrix element
        return real(8. * amp * conj(amp));
      }

      // wave function,
      // ixx true takes anti-particle, false takes particle
      void ixx(std::vector<double> p, int nhel, std::complex<double> fi[6], bool isixx) {
        std::complex<double> chi[2];
        double sqp0p3;
        fi[0] = std::complex<double> (p[0], p[3]);
        fi[1] = std::complex<double> (p[1], p[2]);
        if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0) sqp0p3 = 0.0;
        else sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5);
        if (isixx) chi[0] = std::complex<double> (-sqp0p3, 0.0);
        else chi[0] = std::complex<double> (sqp0p3, 0.0);
        if (sqp0p3 == 0.0) chi[1] = std::complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0);
        else chi[1] = std::complex<double> (nhel * p[1], -p[2])/sqp0p3;
        if (isixx) {
          if (nhel == 1) {
            fi[2] = chi[1];
            fi[3] = chi[0];
            fi[4] = std::complex<double> (0.0, 0.0);
            fi[5] = std::complex<double> (0.0, 0.0);
          } else {
            fi[2] = std::complex<double> (0.0, 0.0);
            fi[3] = std::complex<double> (0.0, 0.0);
            fi[4] = chi[0];
            fi[5] = chi[1];
          }
        } else {
          if(nhel == 1){
            fi[2] = chi[0];
            fi[3] = chi[1];
            fi[4] = std::complex<double> (0.00, 0.00);
            fi[5] = std::complex<double> (0.00, 0.00);
          } else {
            fi[2] = std::complex<double> (0.00, 0.00);
            fi[3] = std::complex<double> (0.00, 0.00);
            fi[4] = chi[1];
            fi[5] = chi[0];
          }
        }
        return;
      }

      // vertices
      void VVS2_3(std::complex<double> V1[], std::complex<double> V2[],
                  std::complex<double> COUP,
                  double M3, double W3, std::complex<double> S3[]) {
        std::complex<double> cI = std::complex<double> (0., 1.);
        std::complex<double> TMP1;
        double P3[4];
        std::complex<double> denom;
        S3[0] = +V1[0] + V2[0];
        S3[1] = +V1[1] + V2[1];
        P3[0] = -S3[0].real();
        P3[1] = -S3[1].real();
        P3[2] = -S3[1].imag();
        P3[3] = -S3[0].imag();
        TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]);
        denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
                      M3 * (M3 - cI * W3));
        S3[2] = denom * cI * TMP1;
      }

      void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
                    std::complex<double> COUP1, std::complex<double> COUP2,
                    double M3, double W3, std::complex<double> V3[]) {

        std::complex<double> cI = std::complex<double> (0., 1.);
        std::complex<double> denom;
        std::complex<double> TMP11;
        double P3[4];
        std::complex<double> TMP14;
        double OM3 = 1./pow(M3, 2);
        V3[0] = +F1[0] + F2[0];
        V3[1] = +F1[1] + F2[1];
        P3[0] = -V3[0].real();
        P3[1] = -V3[1].real();
        P3[2] = -V3[1].imag();
        P3[3] = -V3[0].imag();
        TMP14 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
                 F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
        TMP11 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
                 F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
        denom = 1. / (pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
                      M3 * (M3 - cI * W3));
        V3[2] = COUP2 * denom * - 2. * cI * (OM3 * - 1./2. * P3[0] * (TMP11 + 2. * (TMP14)) +
                      (+1./2. * (F2[4] * F1[2] + F2[5] * F1[3]) + F2[2] * F1[4] + F2[3] * F1[5]));
        V3[3] = COUP2 * denom * - 2. * cI * (OM3 * - 1./2. * P3[1] * (TMP11 + 2. * (TMP14)) +
                      (-1./2. * (F2[5] * F1[2] + F2[4] * F1[3]) + F2[3] * F1[4] + F2[2] * F1[5]));
        V3[4] = COUP2 * denom * 2. * cI * (OM3 * 1./2. * P3[2] * (TMP11 + 2. * (TMP14)) +
                      (+1./2. * cI * (F2[5] * F1[2]) - 1./2. * cI * (F2[4] * F1[3]) - cI *
                       (F2[3] * F1[4]) + cI * (F2[2] * F1[5])));
        V3[5] = COUP2 * denom * 2. * cI * (OM3 * 1./2. * P3[3] * (TMP11 + 2. * (TMP14)) +
                      (+1./2. * (F2[4] * F1[2]) - 1./2. * (F2[5] * F1[3]) - F2[2] * F1[4] + F2[3] * F1[5]));
        V3[2] += COUP1 * denom * - cI * (F2[4] * F1[2] + F2[5] * F1[3] - P3[0] * OM3 * TMP11);
        V3[3] += COUP1 * denom * - cI * (-F2[5] * F1[2] - F2[4] * F1[3] - P3[1] * OM3 * TMP11);
        V3[4] += COUP1 * denom * - cI * (-cI * (F2[5] * F1[2]) + cI * (F2[4] * F1[3]) - P3[2] * OM3 * TMP11);
        V3[5] += COUP1 * denom * - cI * (F2[5] * F1[3] - F2[4] * F1[2] - P3[3] * OM3 * TMP11);
      }

      std::complex<double> VVS3_0(double mass,
                                  std::complex<double> S3[],
                                  std::complex<double> COUP) {
        std::complex<double> TMP15;
        TMP15 = std::complex<double> (0., pow(mass, 2) / 2.);
        return COUP * S3[2] * (TMP15);
      }

      static const int nwavefuncs = 9;
      std::complex<double> w[nwavefuncs][18];

      // Pointer to the model parameters
      Parameters_heft pars;

      // vector with momenta (to be changed each event)
      std::vector < std::vector<double> > pout;
    };

    CPPProcess_P0_Sigma_heft_pp_H_ZZ_4l_heft_gg_epemmupmum MGME;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2020_I1790439);
}
