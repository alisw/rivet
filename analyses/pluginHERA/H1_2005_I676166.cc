// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISLepton.hh"
#include "Rivet/Projections/DISFinalState.hh"


namespace Rivet {


  /// @brief Measurement of beauty production at HERA using events with muons and jets
  class H1_2005_I676166 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2005_I676166);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(Cuts::abspid == PID::MUON), "muons");
      declare(DISKinematics(), "Kinematics");
      const DISLepton dl;
      //declare(dl, "Lepton");
      declare(dl.remainingFinalState(), "FS");

      // DIS events: final state particles boosted to Breit frame then clustered
      // using FastJet KT algorithm with jet radius parameter 1 , pT recombination scheme
      const DISFinalState DISfs(DISFinalState::BoostFrame::BREIT);
      declare(FastJets(DISfs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 1.0,
                      JetAlg::Muons::ALL, JetAlg::Invisibles::NONE, nullptr), "DISjets");

      // Photoproduction events: final state particles in lab frame then clustered
      // using FastJet KT algorithm with jet radius parameter 1 , pT recombination scheme
      const DISFinalState PHOfs(DISFinalState::BoostFrame::LAB);
      declare(FastJets(PHOfs, fastjet::JetAlgorithm::kt_algorithm, fastjet::RecombinationScheme::pt_scheme, 1.0,
                       JetAlg::Muons::ALL, JetAlg::Invisibles::NONE, nullptr), "PHOjets");

      // Photoproduction (Table 4)
      book(_h["PHO_eta_mu"], 1, 1, 1);
      book(_h["PHO_pT_mu"], 2, 1, 1);
      book(_h["PHO_pT_jet"], 3, 1, 1);
      book(_h["PHO_x_obs_gamma"], 4, 1, 1);

      // Electroproduction (Table 5)
      book(_h["DIS_Q^2"], 5, 1, 1);
      book(_h["DIS_x"], 6, 1, 1);
      book(_h["DIS_eta_mu"], 7, 1, 1);
      book(_h["DIS_pT_mu"], 8, 1, 1);
      book(_h["DIS_pT_jet_Breit"], 9, 1, 1);

      isDIS  = false;
      nPHO   = 0;
      nDIS   = 0;
      nVeto0 = 0;
      nVeto1 = 0;
      nVeto2 = 0;
      nVeto3 = 0;
      nVeto4 = 0;
      nVeto5 = 0;
      nVeto6 = 0;
      nVeto7 = 0;

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      isDIS  = false;

      //First Lorentz invariant quantities in Lab frame
      DISKinematics dis = apply<DISKinematics>(event, "Kinematics");
      //const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      double Q2 = dis.Q2()/GeV2;
      double x = dis.x();
      double y = dis.y();

      //after x_g
      //calculated from mass of

      // Separate into DIS and PHO regimes else veto
      if (Q2 < 1 && inRange(y, 0.2, 0.8)) {
        isDIS = false;
        ++nPHO;
      } else if (inRange(Q2, 2.0, 100) && inRange(y, 0.1, 0.7)) {
        isDIS = true;
        ++nDIS;
      } else {
        vetoEvent;
      }
      ++nVeto0;

      // Extract the particles other than the lepton
      Particles particles = apply<FinalState>(event, "FS").particles();

      // Get Lorentz transforms for Breit Boost and Lab Boost
      const LorentzTransform breitboost = dis.boostBreit();
      const LorentzTransform hcmboost = dis.boostHCM(); // Hadron cm system
      const LorentzTransform labboost = breitboost.inverse();


      // If DIS regime, cluster in Breit frame. If PHO regime, cluster in lab frame
      // Apply jet cuts for relevant regime
      Jets DISjets, PHOjets;
      Jets cutJetsDIS, cutJetsPHO;
      Jet pjet1, pjet2;

      if (isDIS) {

        DISjets = apply<FastJets>(event, "DISjets").jetsByPt(Cuts::pT > 6*GeV);
        // Cut on Pseudorapidity in lab frame
        for (vector<int>::size_type i=0; i<DISjets.size(); i++){
           const FourMomentum LabMom = labboost.transform(DISjets[i]);
           // eta cut is in lab frame
           double etaDISJet =abs(LabMom.eta());
           if(etaDISJet < 2.5 ){
              cutJetsDIS.push_back(DISjets[i]);
           }
        }

        // Apply number of jet cuts: in DIS at least one jet is needed
        if (cutJetsDIS.size()<1)		vetoEvent;
        ++nVeto1;

      } else {

        PHOjets = apply<FastJets>(event, "PHOjets").jetsByPt(Cuts::pT > 6*GeV);

        // Cut on Pseudorapidity
        for(vector<int>::size_type i=0; i<PHOjets.size(); i++){
          double etaPHOJet = PHOjets[i].eta();
          if(abs(etaPHOJet) < 2.5 ){
            cutJetsPHO.push_back(PHOjets[i]);
          }
        }

        // Apply number of jet cuts: in phtotproduction at least 2 jets are needed

        if(cutJetsPHO.size()<2)		vetoEvent;
        ++nVeto2;

        // Ensure 1st (2nd) hardest jets have pT > 7 (6) GeV
        pjet1 = cutJetsPHO[0];
        pjet2 = cutJetsPHO[1];

        if(!(pjet1.pT()>7 || pjet2.pT()>7) ){
          vetoEvent;
        }

        ++nVeto3;

      }

      // Apply muon cuts in relevant regime
      const Particles& all_muons = apply<FinalState>(event, "muons").particlesByPt();
      Particles DIS_muons_cut;
      Particles PHO_muons_cut;

      // Apply muon eta / pT cuts
      const Particles DIS_muons = filter_select(all_muons, [](const Particle& m) {
            return m.eta() > -0.75 && m.eta() < 1.15 && m.pT() > 2.5*GeV; });
      const Particles PHO_muons = filter_select(all_muons, [](const Particle& m) {
            return m.eta() > -0.55 && m.eta() < 1.10 && m.pT() > 2.5*GeV; });

      if (isDIS) {
        // Veto event events with no muons
        if (DIS_muons.size() == 0) vetoEvent;
        ++nVeto4;
      }
      else {
        // Veto event events with no muons
        if (PHO_muons.size() == 0) vetoEvent;
        ++nVeto5;
      }

      // Calculate fraction of photon energy entering the hard interaction observable
      double sumJet1 = 0;
      double sumJet2 = 0;
      double sumAllHad = 0;
      double x_gamma;

      if (!isDIS) {
        //const Particles& hadfs = apply<HadronicFinalState>(event, "HFS").particles();
        //for (const Particle& h : hadfs) {
        for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
          Particle& h = particles[ip1];
          FourMomentum hcmMom = hcmboost.transform(h.momentum());
          // Need to change sign: by default hcmboost has gamma* in +z dir,
          // and p in -z dir, but here we need: gamma* -in -z and proton in +z.
          sumAllHad += (hcmMom.E() + hcmMom.pz());
          if ( deltaR(h, pjet1) <= 1.0 ) {
            sumJet1 += (hcmMom.E() + hcmMom.pz());
          }
          else if ( deltaR(h, pjet2) <= 1.0 ) {
            sumJet2 += (hcmMom.E() + hcmMom.pz());
          }
        }

        //sumJet1 = hcmboost.transform(pjet1).E() + hcmboost.transform(pjet1).pz();
        //sumJet2 = hcmboost.transform(pjet2).E() + hcmboost.transform(pjet2).pz();

        x_gamma = (sumJet1 + sumJet2)/sumAllHad;

      }


      // Fill histos
      if (isDIS ) {

        for (const Particle& m : DIS_muons) {
          _h["DIS_pT_mu"] -> fill(m.pT());
          _h["DIS_eta_mu"] -> fill(m.eta());
        }

        _h["DIS_pT_jet_Breit"] -> fill(cutJetsDIS[0].pT()/GeV);
        _h["DIS_Q^2"] -> fill(Q2);
        _h["DIS_x"] -> fill(log10(x));

      } else {

        for (const Particle& m : PHO_muons) {
          _h["PHO_pT_mu"] -> fill(m.pT());
          _h["PHO_eta_mu"] -> fill(m.eta());
        }

        _h["PHO_pT_jet"] -> fill(cutJetsPHO[0].pT());
        _h["PHO_x_obs_gamma"] -> fill(x_gamma);

      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      double normpb = crossSection()/picobarn/sumW();

      scale(_h["PHO_eta_mu"], normpb);
      scale(_h["PHO_pT_mu"], normpb);
      scale(_h["PHO_pT_jet"], normpb);
      scale(_h["PHO_x_obs_gamma"], normpb);

      scale(_h["DIS_Q^2"], normpb);
      scale(_h["DIS_x"], normpb);
      scale(_h["DIS_eta_mu"], normpb);
      scale(_h["DIS_pT_mu"], normpb);
      scale(_h["DIS_pT_jet_Breit"], normpb);

      MSG_DEBUG("Events passing Q2/y cuts   = " << nVeto0 );
      MSG_DEBUG("PHO events   = " << nPHO );
      MSG_DEBUG("DIS events = " << nDIS );
      MSG_DEBUG("DIS Events passing number of jets cuts= " << nVeto1  );
      MSG_DEBUG("PHO Events passing number of jets cuts= " << nVeto2   );
      MSG_DEBUG("PHO Events passing pT jet cuts= " << nVeto3   );
      MSG_DEBUG("DIS Events passing one muon cut= " << nVeto4   );
      MSG_DEBUG("PHO Events passing one muon cut= " << nVeto5   );
      MSG_DEBUG("DIS Events passing muon in jet cut  = " << nVeto6 );
      MSG_DEBUG("PHO Events passing muon in jet cut  = " << nVeto7 );

    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, Profile1DPtr> _p;
    map<string, CounterPtr> _c;
    /// @}

    bool isDIS;
    int nPHO, nDIS;
    int  nVeto0, nVeto1, nVeto2, nVeto3, nVeto4, nVeto5, nVeto6, nVeto7;

  };


  RIVET_DECLARE_PLUGIN(H1_2005_I676166);

}
