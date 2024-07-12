// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"

namespace Rivet {

  /// @brief ATLAS ZZ --> 4l and 2l2v
  class ATLAS_2016_I1494075 : public Analysis {
  public:
    /// Constructor
    //
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_I1494075);
    /// @name Analysis methods
    /// @{
    /// Book histograms and initialise projections before the run
    void init() {
        _mode = 0;
        if (getOption("LMODE") == "2L2NU") _mode = 2;
        if (getOption("LMODE") == "4L") _mode = 1;

        PromptFinalState prompt_photons(Cuts::abspid == PID::PHOTON);
        PromptFinalState prompt_ele(Cuts::abspid == PID::ELECTRON);
        PromptFinalState prompt_mu(Cuts::abspid == PID::MUON);

        // Wide lepton cuts which cover both channels and are used for the jet veto.
        Cut dressedele_cuts = (Cuts::abseta < 4.9) && (Cuts::pT > 7*GeV);
        Cut dressedmu_cuts = (Cuts::abseta < 2.7) && (Cuts::pT > 7*GeV);
        const DressedLeptons dressedelectrons(prompt_photons, prompt_ele, 0.1, dressedele_cuts);
        const DressedLeptons dressedmuons(prompt_photons, prompt_mu, 0.1, dressedmu_cuts);

        declare(dressedelectrons, "electrons");
        declare(dressedmuons, "muons");

        VisibleFinalState vfs;
        VetoedFinalState jetinput(vfs);
        jetinput.addVetoOnThisFinalState(dressedmuons);

        if (_mode != 1)  declare(InvisibleFinalState(true), "MET");

        FastJets fastjets(jetinput, FastJets::ANTIKT, 0.4);
        declare (fastjets, "Jets");

        // ZZ to four leptons channel
        book(_h["leading_ll_pt"], 2, 1, 1);
        book(_h["Njets"], 3, 1, 1);
        book(_h["leading_ll_phi"], 4, 1, 1);
        book(_h["ZZ_rapidity"], 5, 1, 1);
        // ZZ to lvlv channel
        book(_h["dilepton_pt"], 6, 1, 1);
        book(_h["llphi_lvchannel"], 7, 1, 1);
        book(_h["mzz_lvchannel"], 8, 1, 1);

    }

    struct Zstate : public ParticlePair {
      Zstate() { }
      Zstate(ParticlePair _particlepair) : ParticlePair(_particlepair) { }
      FourMomentum mom() const { return first.momentum() + second.momentum();}
      double rapid() const { return ((first.momentum() + second.momentum()).rapidity()); }
      double dphi() const { return deltaPhi(first.ptvec(), second.ptvec()); }
    };


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // find out how many good jets
      //Jets jets = apply<FastJets>(event, "Jets").jetsByPt();
      Particles cand_e = apply<DressedLeptons>(event, "electrons").particlesByPt();
      Particles cand_mu = apply<DressedLeptons>(event, "muons").particlesByPt();

      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < 4.5 && Cuts::pT>25*GeV);
      idiscardIfAnyDeltaRLess(jets, cand_e, 0.3);

      //llvv channel
      if (_mode != 1) {

        // jet veto
        if (jets.empty()){

          Particles selected_pair;
          Vector3 met_vec;

          const FinalState& metfs = apply<InvisibleFinalState>(event, "MET");
          for (const Particle& p : metfs.particles())  met_vec += p.mom().perpVec();

          for ( const Particle& mu : cand_mu ) {
            if (mu.pT() > 25*GeV && mu.abseta() < 2.5 ) {
              selected_pair.push_back(mu);
            }
          }
          for ( const Particle& e : cand_e ) {
            if (e.pT() > 25*GeV && e.abseta() < 2.5 ) {
              selected_pair.push_back(e);
            }
          }

          // selections on pair.
          if ((selected_pair.size() == 2) // exactly two leptons
              && (selected_pair[0].abspid() == selected_pair[1].abspid()) //same flavour
          && (deltaR(selected_pair[0],selected_pair[1]) > 0.3) // deltaR > 0.3
          && (selected_pair[0].pid() * selected_pair[1].pid() < 0)) {  // opposite sign

            //ensure the first one is leading lepton
            if (selected_pair[0].momentum().pT() < selected_pair[1].momentum().pT()) {
              std::swap(selected_pair[0], selected_pair[1]);
            }

            //ZZ four momentum, three momentum if from invisble final states
            const FourMomentum Z_1_mom = selected_pair[0].momentum() + selected_pair[1].momentum();

            const double axial_Etmiss = -1.0*met_vec.mod()*cos(deltaPhi(met_vec, Z_1_mom.ptvec()));

            double pT_balance = fabs( (met_vec.mod() - Z_1_mom.pT()) /Z_1_mom.pT() );
            if (axial_Etmiss > 90*GeV && pT_balance < 0.4 && inRange(Z_1_mom.mass(), 76*GeV, 106*GeV) ) {

              double mz_pdg2 = 91.1876*91.1876*GeV*GeV;
              // transverse mass
              double mTrans = sqrt(
                 sqr((sqrt(Z_1_mom.pT()*Z_1_mom.pT() + mz_pdg2) + sqrt(met_vec.mod2() + mz_pdg2)))
                 - (Z_1_mom.ptvec() + met_vec.perpVec()).mod2()
              );
              _h["mzz_lvchannel"]->fill(mTrans/GeV);
              _h["llphi_lvchannel"]->fill(deltaPhi(selected_pair[0].momentum().ptvec(), selected_pair[1].momentum().ptvec()));
              _h["dilepton_pt"]->fill(Z_1_mom.pT()/GeV);

            }
          }
        }
      }

      //for llll
      if (_mode != 2) {

        ///////////
        // Insert selected muons then electrons into the lepton 4l final state
        ///////////
        vector<DressedLepton> leptonsFS_sel4l;
        leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), cand_mu.begin(), cand_mu.end() );
        leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), cand_e.begin(), cand_e.end() );

        ////////////
        // Cut dR>0.2 between all leptons
        Particles n_parts;
        for (const DressedLepton& l1 : leptonsFS_sel4l) {
          bool isolated = true;
          for (DressedLepton& l2 : leptonsFS_sel4l){
                  const double fourL_dR = deltaR(l1, l2);
                  if (fourL_dR < 0.2 && !isSame(l1, l2)) {
              isolated = false;
              break;
            }
          }
          if (isolated) n_parts.push_back(l1);
        }

        double totalCharge = 0;
        for (const Particle& p : n_parts) totalCharge += p.pid();

        if (n_parts.size() == 4 && totalCharge == 0 ) {

          Zstate lead_Z, sub_Z;
          identifyZstates(lead_Z, sub_Z, n_parts);
          if (lead_Z.mom().pT() < sub_Z.mom().pT()) {
            std::swap(lead_Z, sub_Z);
          }

          vector<DressedLepton> lepton4l;
          lepton4l.insert( lepton4l.end(), leptonsFS_sel4l.begin(), leptonsFS_sel4l.end() );
          std::sort(lepton4l.begin(), lepton4l.end(), [](const DressedLepton& l1, const DressedLepton& l2) {
            return (l1.abseta() > l2.abseta());
          });
          if (lead_Z.first.abspid() == 11 && sub_Z.first.abspid() == 11) {
            if (lepton4l[1].abseta() > 2.5 || lepton4l[2].abseta() > 2.5 || lepton4l[3].abseta() > 2.5) vetoEvent;
          }
          else if (lead_Z.first.abspid() != sub_Z.first.abspid()) {
            if (lead_Z.first.abspid() == 11) {
               if (std::min(lead_Z.first.abseta(), lead_Z.second.abseta()) > 2.5)  vetoEvent;
            }
            if (lead_Z.first.abspid() == 13) {
               if (std::min(sub_Z.first.abseta(), sub_Z.second.abseta()) > 2.5)  vetoEvent;
            }
          }

          double m_Z1      = lead_Z.mom().mass();
          double m_Z2      = sub_Z.mom().mass();
          double lead_Z_Pt = lead_Z.mom().pT();
          double lead_dPhi  = lead_Z.dphi();

          double ZZ_rap    = fabs(lead_Z.rapid() - sub_Z.rapid());

          //Z mass selections
          if ( inRange(m_Z1, 66*GeV, 116*GeV) && inRange(m_Z2, 66*GeV, 116*GeV) ) {
            _h["leading_ll_pt"]->fill(lead_Z_Pt/GeV);
            _h["leading_ll_phi"]->fill(lead_dPhi);
            _h["ZZ_rapidity"]->fill(ZZ_rap);
            _h["Njets"]->fill(jets.size());
          }
        }
      }
    };



    /// Normalise histograms etc., after the run
    void finalize() {
      // histo1D is divided by bin width when converting to scatter2D so no need of further normalisation for this one.
      const double sf  = crossSectionPerEvent()/femtobarn;
      const double sf2 = crossSectionPerEvent()/picobarn;
      // 4l is divided by branching ratio to make a ZZ cross section
      const double br = (3.3632 + 3.3662)/100.;
      scale(_h["leading_ll_pt"],  sf/sqr(br));
      scale(_h["Njets"],          sf/sqr(br));
      scale(_h["leading_ll_phi"], sf2/sqr(br));
      scale(_h["ZZ_rapidity"],    sf2/sqr(br));
      // llvv is cross section for a single flavour. The analysis assumes we have run on e and mu,
      // so divides by 2.
      scale(_h["dilepton_pt"], sf/2.0);
      scale(_h["mzz_lvchannel"], sf/2.0);
      scale(_h["llphi_lvchannel"],sf/2.0);

    }

  private:

    void identifyZstates(Zstate& Z1, Zstate& Z2, const Particles& n_parts);
    const double Zmass = 91.1876*GeV; // GeV
    map<string, Histo1DPtr> _h;
    size_t _mode;

  };


  void ATLAS_2016_I1494075::identifyZstates(Zstate& Z1, Zstate& Z2, const Particles& n_parts){


    // first find the lepton types
    Particles part_pos_el, part_neg_el, part_pos_mu, part_neg_mu;
    for (const Particle& l : n_parts) {
      if (l.abspid() == PID::ELECTRON) {
        if (l.pid() < 0) part_neg_el.push_back(l);
        if (l.pid() > 0) part_pos_el.push_back(l);
      }
      else if (l.abspid() == PID::MUON) {
        if (l.pid() < 0) part_neg_mu.push_back(l);
        if (l.pid() > 0) part_pos_mu.push_back(l);
      }
    }

    //4e/4mu channel, pairing ambiguity
    if (part_neg_el.size() == 2 || part_neg_mu.size() == 2) {
      Zstate Zcand1, Zcand2, Zcand3, Zcand4;
      if (part_neg_el.size() == 2) {
        Zcand1 = Zstate( ParticlePair( part_neg_el[0], part_pos_el[0] ) );
        Zcand2 = Zstate( ParticlePair( part_neg_el[0], part_pos_el[1] ) );
        Zcand3 = Zstate( ParticlePair( part_neg_el[1], part_pos_el[0] ) );
        Zcand4 = Zstate( ParticlePair( part_neg_el[1], part_pos_el[1] ) );
      } else {
        Zcand1 = Zstate( ParticlePair( part_neg_mu[0], part_pos_mu[0] ) );
        Zcand2 = Zstate( ParticlePair( part_neg_mu[0], part_pos_mu[1] ) );
        Zcand3 = Zstate( ParticlePair( part_neg_mu[1], part_pos_mu[0] ) );
        Zcand4 = Zstate( ParticlePair( part_neg_mu[1], part_pos_mu[1] ) );
      }

      // pairing should be |1 + 4| and |2 + 3| in mass order
      double V1, V2;
      V1 = fabs( Zcand1.mom().mass() - Zmass ) + fabs( Zcand4.mom().mass() - Zmass);
      V2 = fabs( Zcand2.mom().mass() - Zmass ) + fabs( Zcand3.mom().mass() - Zmass);

      if (V1 > V2) {
        Z1 = Zcand2;
        Z2 = Zcand3;
      }
      else {
        Z1 = Zcand1;
        Z2 = Zcand4;
      }
      //2e2mu
    }
    else if (part_neg_el.size() == 1 && part_neg_mu.size() == 1) {
      Z1 = Zstate( ParticlePair( part_neg_mu[0], part_pos_mu[0] ) );
      Z2 = Zstate( ParticlePair( part_neg_el[0], part_pos_el[0] ) );
    }

  }

  RIVET_DECLARE_PLUGIN(ATLAS_2016_I1494075);
}
