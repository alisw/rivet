// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/IndirectFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/Smearing.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// Recursive jigsaw chargino-neutralino search with 2 or 3 charged leptons in 36/fb of 13 TeV pp
  ///
  /// @author Derek Yeung, Andy Buckley
  class ATLAS_2018_I1676551 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2018_I1676551);


    /// Analysis initialization
    void init() {

      PromptFinalState electrons(Cuts::abspid == PID::ELECTRON);
      SmearedParticles recoelectrons(electrons, ELECTRON_EFF_ATLAS_RUN2_MEDIUM, ELECTRON_SMEAR_ATLAS_RUN2);
      declare(recoelectrons, "Electrons");

      PromptFinalState muons(Cuts::abspid == PID::MUON);
      SmearedParticles recomuons(muons, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2);
      declare(recomuons, "Muons");

      FastJets jets4(IndirectFinalState(Cuts::open()), FastJets::ANTIKT, 0.4);
      SmearedJets recojets4(jets4, JET_SMEAR_CMS_RUN2, JET_BTAG_EFFS(0.77, 1/6., 1/134.));
      declare(recojets4, "Jets");

      MissingMomentum met(FinalState(Cuts::abseta < 4.9));
      SmearedMET recomet(met, MET_SMEAR_ATLAS_RUN2);
      declare(recomet, "MET");

      // Cutflow Setup for 2l-High
      const strings cfnames1 = {"Trigger matching & 2 signal leptons", "Preselection",
                                "Fraction 1 > 0.8", "Fraction 2 < 0.05",
                                "Delta Phi in [0.3, 2.9]", "H_PP_4,1 > 800 GeV"};
      _cutflow2l[0].addCutflow("ATLAS_2018_I1676551 SR EW", cfnames1);

      // Cutflow Setup for 2l-Int
      const strings cfnames2 = {"Trigger matching & 2 signal leptons", "Preselection",
                                "Fraction 1 > 0.8", "Fraction 2 < 0.05",
                                "Delta Phi in [0.3, 2.6]", "H_PP_4,1 > 600 GeV"};
      _cutflow2l[1].addCutflow("ATLAS_2018_I1676551 SR EW", cfnames2);

      // Cutflow Setup for 2l-Low
      const strings cfnames3 = {"Trigger matching & 2 signal leptons", "Preselection",
                                "Fraction 1 in [0.35, 0.6]", "Fraction 2 < 0.05",
                                "Min Delta Phi > 2.4", "H_PP_4,1 > 400 GeV"};
      _cutflow2l[2].addCutflow("ATLAS_2018_I1676551 SR EW", cfnames3);

      // Cutflow Setup for 2l-ISR
      const strings cfnames4 = {"Trigger matching & 2 signal leptons", "Preselection",
                                "m_Z in [80, 100] GeV", "m_J in [50, 110] GeV",
                                "Delta Phi > 2.8", "R_ISR in [0.4, 0.75]", "p_CM_T_ISR > 180 GeV",
                                "p_CM_T_I > 100 GeV", "p_CM_T < 30 GeV"};
      _cutflow2lISR.addCutflow("ATLAS_2018_I1676551 SR EW", cfnames4);

      // Cutflow Setup for 3l-High
      const strings cfnames5 = {"Trigger matching & 3 signal leptons", "Preselection",
                                "m_ll in [75,105] GeV", "m_W_T > 150 GeV",
                                "Fraction 1 > 0.75", "Fraction 2 < 0.8",
                                "H_PP_31 > 500 GeV", "Fraction 3 < 0.2"};
      _cutflow3l[0].addCutflow("ATLAS_2018_I1676551 SR EW", cfnames5);

      // Cutflow Setup for 3l-Int
      const strings cfnames6 = {"Trigger matching & 3 signal leptons", "Preselection",
                                "m_ll in [75,105] GeV", "m_W_T > 130 GeV",
                                "Fraction 1 > 0.8", "Fraction 2 < 0.75",
                                "H_PP_31 > 450 GeV", "Fraction 3 < 0.15"};
      _cutflow3l[1].addCutflow("ATLAS_2018_I1676551 SR EW", cfnames6);

      // Cutflow Setup for 3l-Low
      const strings cfnames7 = {"Trigger matching & 3 signal leptons", "Preselection",
                                "m_ll in [75,105] GeV", "m_W_T > 100 GeV",
                                "Fraction 1 > 0.9", "H_PP_31 > 250 GeV", "Fraction 2 < 0.05"};
      _cutflow3l[2].addCutflow("ATLAS_2018_I1676551 SR EW", cfnames7);

      // Cutflow Setup for 3l-ISR
      const strings cfnames8 = {"Trigger matching & 3 signal leptons", "Preselection",
                                "m_ll in [75, 105] GeV", "m_W_T > 100 GeV",
                                "Delta Phi > 2.0", "R_ISR in [0.55, 1.0]", "p_CM_T_ISR > 100 GeV",
                                "p_CM_T_I > 80 GeV", "p_CM_T < 25 GeV"};
      _cutflow3lISR.addCutflow("ATLAS_2018_I1676551 SR EW", cfnames8);
    }


    // Per-event analysis
    void analyze(const Event& event) {

      _cutflow2l[0].fillinit();
      _cutflow2l[1].fillinit();
      _cutflow2l[2].fillinit();
      _cutflow2lISR.fillinit();
      _cutflow3l[0].fillinit();
      _cutflow3l[1].fillinit();
      _cutflow3l[2].fillinit();
      _cutflow3lISR.fillinit();

      // Obtain Electrons, Muons and Jets
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particlesByPt(Cuts::pT > 10*GeV && Cuts::abseta < 2.47);
      Particles muons = apply<ParticleFinder>(event, "Muons").particlesByPt(Cuts::pT > 10*GeV && Cuts::abseta < 2.4);
      Particles leptons = sortByPt(elecs + muons);
      Jets jets = apply<SmearedJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);

      // Discard jets within DR = 0.4 of prompt leptons
      idiscardIfAnyDeltaRLess(jets, leptons, 0.4);

      // 2-lepton High (n=0), Int (n=1) and Low (n=2) Selection
      for (int n=0; n<3; ++n) {

        while (true) {

          // Require the leading electron and the leading muon to have pT > 25GeV
          if (elecs.size() != 0 && elecs[0].pT() < 25*GeV) break;
          if (muons.size() != 0 && muons[0].pT() < 25*GeV) break;

          // Obtain the missing transverse momentum vector
          Vector3 EtMissX = apply<SmearedMET>(event,"MET").vectorPt();
          Vector3 EtMiss  = -EtMissX;

          // Require only 2 opposite charged, same flavoured leptons
          double n_leptons = leptons.size();
          if (n_leptons != 2) break;
          if (leptons[0].abspid() != leptons[1].abspid()) break;
          if (leptons[0].charge() * leptons[1].charge() >= 0) break;

          // Obtain the 4-momenta of leptons and implement requirements on them
          vector<FourMomentum> lepton;
          for (int l = 0; l < 2; ++l) lepton.push_back(leptons[l].mom());

          double p_l1_T = lepton[0].pT();
          if (p_l1_T < 25*GeV) break;

          double p_l2_T = lepton[1].pT();
          if (p_l2_T < 25*GeV) break;

          _cutflow2l[n].fill(1);

          // Requirement on m_ll
          double m_ll = (lepton[0]+lepton[1]).mass();
          if (!inRange(m_ll, 80*GeV, 100*GeV)) break;

          // Obtain the 4-momenta of jets and implement requirements on them
          double n_jets = jets.size();
          if (n != 2) {
            if (n_jets < 2) break;
            if (any(jets,hasBTag())) break;
          } else{
            if (n_jets != 2) break;
            if (any(jets,hasBTag())) break;
          }

          vector<FourMomentum> jet;
          for (int l = 0; l < 2; ++l) jet.push_back(jets[l].momentum());

          double p_j1_T = jet[0].pT();
          if (p_j1_T < 30*GeV) break;

          double p_j2_T = jet[1].pT();
          if (p_j2_T < 30*GeV) break;

          double m_jj = (jet[0]+jet[1]).mass();
          if (n != 2) {
            if (m_jj < 60*GeV || m_jj > 100*GeV) break;
          } else{
            if (m_jj < 70*GeV || m_jj > 90*GeV) break;
          }

          // Pre-selection requirements cut
          _cutflow2l[n].fill(2);

          // Invisible mass JR and Invisible Rapidity JR to obtain Invisible System 4-momentum
          FourMomentum P_V = lepton[0] + lepton[1] + jet[0] + jet[1];
          double M_I = sqrt(P_V.mass2() - 4*m_ll*m_jj);
          double M_V = P_V.mass();
          double p_I_z = P_V.pz() * sqrt(EtMiss.dot(EtMiss) + sqr(M_I)) / sqrt(P_V.pT2() + sqr(M_V));

          FourMomentum P_I;
          P_I.setPM(EtMiss.x(), EtMiss.y(), p_I_z, M_I);

          // Lorentz boost to CM frame
          LorentzTransform LT = LorentzTransform::mkFrameTransform(P_I+P_V);

          FourMomentum P_F_I; //4-momentum of invisible system in CM frame
          FourMomentum P_F_Va; //4-momentum of visible system a in CM frame
          FourMomentum P_F_Vb; //4-momentum of visible system b in CM frame
          FourMomentum P_F_V; //4-momentum of visible system in CM frame
          FourMomentum P_F_Va1; //4-momentum of visible system a1 in CM frame
          FourMomentum P_F_Va2; //4-momentum of visible system a2 in CM frame
          FourMomentum P_F_Vb1; //4-momentum of visible system b1 in CM frame
          FourMomentum P_F_Vb2; //4-momentum of visible system b2 in CM frame
          if ((jet[0]+jet[1]).mass() > (lepton[0]+lepton[1]).mass()) {
            P_F_I = LT.transform(P_I);
            P_F_Va = LT.transform(jet[0]+jet[1]);
            P_F_Vb = LT.transform(lepton[0]+lepton[1]);
            P_F_V = LT.transform(P_V);
            //
            P_F_Va1 = LT.transform(jet[0]);
            P_F_Va2 = LT.transform(jet[1]);
            P_F_Vb1 = LT.transform(lepton[0]);
            P_F_Vb2 = LT.transform(lepton[1]);
          } else {
            P_F_I = LT.transform(P_I);
            P_F_Va = LT.transform(lepton[0]+lepton[1]);
            P_F_Vb = LT.transform(jet[0]+jet[1]);
            P_F_V = LT.transform(P_V);
            //
            P_F_Va1 = LT.transform(lepton[0]);
            P_F_Va2 = LT.transform(lepton[1]);
            P_F_Vb1 = LT.transform(jet[0]);
            P_F_Vb2 = LT.transform(jet[1]);
          }

          // Obtain variables defined in the Contra-boost Invariant JR
          double M2c = 2*((P_F_Va.E())*(P_F_Vb.E())+(P_F_Va.p3()).dot(P_F_Vb.p3()));
          double m2Va = sqr(P_F_Va.mass());
          double m2Vb = sqr(P_F_Vb.mass());
          double ka = m2Va-m2Vb+M2c-2*sqrt(m2Va)*sqrt(m2Vb);
          double kb = m2Vb-m2Va+M2c-2*sqrt(m2Va)*sqrt(m2Vb);
          double kn = ka*m2Va - kb*m2Vb + M2c*(kb-ka)/2 + (0.5)*sqrt(pow(ka+kb,2)*(pow(M2c,2)-4*m2Va*m2Vb));
          double k2d = sqr(ka)*m2Va+sqr(kb)*m2Vb+ka*kb*M2c;
          double k = kn/k2d;
          double ca = (1+k*ka)/2;
          double cb = (1+k*kb)/2;
          double c = 0.5*(P_F_V.E()+sqrt(pow(P_F_V.E(),2)+pow(M_I,2)-pow(P_V.mass(),2)))/(ca*(P_F_Va.E())+cb*(P_F_Vb.E()));

          // Apply Contra-boost Invariant JR to obtain invisible particles' 4-momenta
          double p_F_Iax = (P_F_Va.px())*(c*ca-1)-(P_F_Vb.px())*c*cb;
          double p_F_Iay = (P_F_Va.py())*(c*ca-1)-(P_F_Vb.py())*c*cb;
          double p_F_Iaz = (P_F_Va.pz())*(c*ca-1)-(P_F_Vb.pz())*c*cb;
          double p_F_Ibx = (P_F_Vb.px())*(c*cb-1)-(P_F_Va.px())*c*ca;
          double p_F_Iby = (P_F_Vb.py())*(c*cb-1)-(P_F_Va.py())*c*ca;
          double p_F_Ibz = (P_F_Vb.pz())*(c*cb-1)-(P_F_Va.pz())*c*ca;
          double E_F_Ia = (c*ca-1)*(P_F_Va.E())+c*cb*(P_F_Vb.E());
          double E_F_Ib = (c*cb-1)*(P_F_Vb.E())+c*ca*(P_F_Va.E());

          FourMomentum P_F_Ia;
          FourMomentum P_F_Ib;
          P_F_Ia.setPE(p_F_Iax,p_F_Iay,p_F_Iaz,E_F_Ia);
          P_F_Ib.setPE(p_F_Ibx,p_F_Iby,p_F_Ibz,E_F_Ib);

          // Lorentz boost from the CM frame to the P_a and P_b frame
          LorentzTransform LTPa=LorentzTransform::mkFrameTransform(P_F_Va+P_F_Ia);
          LorentzTransform LTPb=LorentzTransform::mkFrameTransform(P_F_Vb+P_F_Ib);

          // min(H_Pa_11,H_Pb_11)/min(H_Pa_21,H_Pb,21) requirement
          if (n != 2) {
            double H_Pa_11 = (LTPa.transform(P_F_Va1+P_F_Va2)).p()+(LTPa.transform(P_F_Ia)).p();
            double H_Pb_11 = (LTPb.transform(P_F_Vb1+P_F_Vb2)).p()+(LTPb.transform(P_F_Ib)).p();
            double H_Pa_21 = (LTPa.transform(P_F_Va1)).p()+(LTPa.transform(P_F_Va2)).p()+(LTPa.transform(P_F_Ia)).p();
            double H_Pb_21 = (LTPa.transform(P_F_Vb1)).p()+(LTPb.transform(P_F_Vb2)).p()+(LTPb.transform(P_F_Ib)).p();
            vector<double> V1 = {H_Pa_11,H_Pb_11};
            vector<double> V2 = {H_Pa_21,H_Pb_21};
            double fraction1 = min(V1)/min(V2);
            if (fraction1 < 0.8) break;
          } else {
            double H_PP_41 = P_F_Va1.p()+P_F_Va2.p()+P_F_Vb1.p()+P_F_Vb2.p()+(P_F_Ia+P_F_Ib).p();
            double H_PP_11 = (P_F_Va1+P_F_Va2+P_F_Vb1+P_F_Vb2).p()+(P_F_Ia+P_F_Ib).p();
            double fraction1 = H_PP_11/H_PP_41;
            if (fraction1 < 0.35 || fraction1 > 0.6) break;
          }
          _cutflow2l[n].fill(3);

          // Lorentz boost from the CM frame to the lab frame
          LorentzTransform LTR = LorentzTransform::mkObjTransform(P_I+P_V);

          // p_lab_T_PP/(p_lab_T_PP+p_lab_T_41) requirement
          double p_lab_T_PP = (LTR.transform(P_F_Va1+P_F_Va2+P_F_Vb1+P_F_Vb2+P_F_Ia+P_F_Ib)).pT();
          double H_PP_T_41 = P_F_Va1.pT()+P_F_Va2.pT()+P_F_Vb1.pT()+P_F_Vb2.pT()+(P_F_Ia+P_F_Ib).pT();
          double fraction2 = p_lab_T_PP/(p_lab_T_PP+H_PP_T_41);
          if (fraction2 > 0.05) break;
          _cutflow2l[n].fill(4);

          // Delta-Phi requirement
          if (n == 0) {
            Vector3 vector1 = (P_F_Va+P_F_Ia).betaVec();
            Vector3 vector2 = (LTPa.transform(P_F_Va)).betaVec();
            if (deltaPhi(vector1,vector2) < 0.3 || deltaPhi(vector1,vector2) > 2.8) break;
            Vector3 vector3 = (P_F_Vb+P_F_Ib).betaVec();
            Vector3 vector4 = (LTPb.transform(P_F_Vb)).betaVec();
            if (deltaPhi(vector3,vector4) < 0.3 || deltaPhi(vector3,vector4) > 2.8) break;
          } else if (n == 1) {
            Vector3 vector1 = (P_F_Va+P_F_Ia).betaVec();
            Vector3 vector2 = (LTPa.transform(P_F_Va)).betaVec();
            if (deltaPhi(vector1,vector2) < 0.6 || deltaPhi(vector1,vector2) > 2.6) break;
            Vector3 vector3 = (P_F_Vb+P_F_Ib).betaVec();
            Vector3 vector4 = (LTPb.transform(P_F_Vb)).betaVec();
            if (deltaPhi(vector3,vector4) < 0.6 || deltaPhi(vector3,vector4) > 2.6) break;
          } else {
            Vector3 j1 = jet[0].p3();
            Vector3 j2 = jet[1].p3();
            double delta1 = deltaPhi(j1,EtMiss);
            double delta2 = deltaPhi(j2,EtMiss);
            double delta = min(delta1,delta2);
            if (delta < 2.4) break;
          }
          _cutflow2l[n].fill(5);

          // H_PP_41 requirement
          double H_PP_41 = P_F_Va1.p()+P_F_Va2.p()+P_F_Vb1.p()+P_F_Vb2.p()+(P_F_Ia+P_F_Ib).p();
          if (n==0) {
            if (H_PP_41 < 800*GeV) break;
          } else if (n==1) {
            if (H_PP_41 < 600*GeV) break;
          } else {
            if (H_PP_41 < 400*GeV) break;
          }
          _cutflow2l[n].fill(6);

          break;
        }
      }


      // 2-lepton ISR Selection
      while (true) {

        // Require the leading electron and the leading muon to have pT > 25 GeV
        if (elecs.size() != 0 && elecs[0].pT() < 25*GeV) break;
        if (muons.size() != 0 && muons[0].pT() < 25*GeV) break;

        // Require only 2 opposite charged, same flavoured leptons
        double n_leptons = leptons.size();
        if (n_leptons != 2) break;
        if (leptons[0].abspid() != leptons[1].abspid()) break;
        if (leptons[0].charge()*leptons[1].charge() >= 0) break;

        // Obtain the 4-momenta of leptons and implement requirements on them
        vector<FourMomentum> lepton;
        double M;
        for (int n = 0; n < 2; ++n) {
          lepton.push_back(leptons[n].mom());
          M = lepton[n].mass();
          lepton[n].setPz(0);
          lepton[n].setE(sqrt(pow(lepton[n].px(),2)+pow(lepton[n].py(),2)+pow(M,2)));
        }
        if (lepton[0].pT() < 25*GeV || lepton[1].pT() < 25*GeV) break;

        _cutflow2lISR.fill(1);

        // Obtain the 4-momenta of leptons and implement requirements on them
        double n_jets = jets.size();
        if (n_jets < 3 || n_jets > 4) break;
        if (any(jets, hasBTag())) break;

        vector<FourMomentum> jet;
        for (int n = 0; n < n_jets; ++n) {
          jet.push_back(jets[n].momentum());
          M = jet[n].mass();
          jet[n].setPz(0);
          jet[n].setE(sqrt(pow(jet[n].px(),2)+pow(jet[n].py(),2)+pow(M,2)));
        }
        if (jet[0].pT() < 30*GeV || jet[1].pT() < 30*GeV) break;

        Vector3 EtMissX = apply<SmearedMET>(event,"MET").vectorPt();
        Vector3 EtMiss = -EtMissX;

        // Obtain the missing transverse momentum vector
        FourMomentum EtMiss1;
        EtMiss1.setPM(EtMiss.x(),EtMiss.y(),EtMiss.z(),sqrt(EtMiss.dot(EtMiss)));

        // Combinatoric Minimization JR
        vector<double> f;
        int indexToReturn = 0;
        int indexValue = 0;
        int newValue = 0;
        if (n_jets == 3) {
          f.push_back((EtMiss1+jet[1]+jet[2]).p());
          f.push_back((EtMiss1+jet[0]+jet[2]).p());
          f.push_back((EtMiss1+jet[0]+jet[1]).p());
          f.push_back((EtMiss1+jet[0]).p());
          f.push_back((EtMiss1+jet[1]).p());
          f.push_back((EtMiss1+jet[2]).p());

          for (int i = 0; i < 6; i++) {
            newValue = f[i];
            if (newValue >= indexValue) {
              indexToReturn = i;
              indexValue = newValue;
            }
          }

          if (indexToReturn == 3 || indexToReturn == 4 || indexToReturn == 5) break;
        }
        else {
          f.push_back((EtMiss1+jet[0]+jet[1]).p());
          f.push_back((EtMiss1+jet[0]+jet[2]).p());
          f.push_back((EtMiss1+jet[0]+jet[3]).p());
          f.push_back((EtMiss1+jet[1]+jet[2]).p());
          f.push_back((EtMiss1+jet[1]+jet[3]).p());
          f.push_back((EtMiss1+jet[2]+jet[3]).p());
          f.push_back((EtMiss1+jet[0]).p());
          f.push_back((EtMiss1+jet[1]).p());
          f.push_back((EtMiss1+jet[2]).p());
          f.push_back((EtMiss1+jet[3]).p());
          f.push_back((EtMiss1+jet[1]+jet[2]+jet[3]).p());
          f.push_back((EtMiss1+jet[0]+jet[2]+jet[3]).p());
          f.push_back((EtMiss1+jet[0]+jet[1]+jet[3]).p());
          f.push_back((EtMiss1+jet[0]+jet[1]+jet[2]).p());

          for (int i = 0; i < 14; i++) {
            newValue = f[i];
            if (newValue >= indexValue) {
              indexToReturn = i;
              indexValue = newValue;
            }
          }

          if (indexToReturn == 6 || indexToReturn == 7 || indexToReturn == 8 || indexToReturn == 9 ||
              indexToReturn == 10 || indexToReturn == 11 || indexToReturn == 12 || indexToReturn == 13) break;
        }

        _cutflow2lISR.fill(2);

        // g1 and g2 are the set of jets belonging to the ISR and signal system respectively
        vector<vector<double>> g1, g2;
        if (n_jets == 3) {
          g2 = {{1,2},{0,2},{0,1}};
        } else {
          g1 = {{2,3},{1,3},{1,2},{0,3},{0,2},{0,1}};
          g2 = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
        }

        // m_Z requirement
        double m_Z = ((leptons[0]).mom()+(leptons[1]).mom()).mass();
        if (m_Z < 80*GeV || m_Z > 100*GeV) break;
        _cutflow2lISR.fill(3);

        // m_J requirement
        double m_J = ((jets[g2[indexToReturn][0]]).momentum()+(jets[g2[indexToReturn][1]]).momentum()).mass();
        if (m_J < 50*GeV || m_J > 110*GeV) break;
        _cutflow2lISR.fill(4);

        // Compute variables delta_phi, R_ISR, P_T_ISR, p_T_I, p_T and implement their requirements
        double p_T_ISR;
        double p_T_I;
        double p_T;
        double R_ISR;
        double delta_phi;
        if (n_jets == 3) {
          FourMomentum CM = EtMiss1+jet[indexToReturn]+jet[g2[indexToReturn][0]]+jet[g2[indexToReturn][1]]+lepton[0]+lepton[1];
          LorentzTransform LT=LorentzTransform::mkFrameTransform(CM);
          p_T_ISR = (LT.transform(jet[indexToReturn])).pT();
          p_T_I = (LT.transform(EtMiss1)).pT();
          p_T = (CM).pT();

          Vector3 p_I = (LT.transform(EtMiss1)).p3();
          Vector3 p_S = (LT.transform(EtMiss1+jet[g2[indexToReturn][0]]+jet[g2[indexToReturn][1]]+lepton[0]+lepton[1])).p3();
          R_ISR = (p_S.dot(p_I))/(p_S.dot(p_S));
          delta_phi = deltaPhi((LT.transform(EtMiss1)).p3(),(LT.transform(jet[indexToReturn])).p3());

        } else {

          FourMomentum CM = EtMiss1+jet[g1[indexToReturn][0]]+jet[g1[indexToReturn][1]]+jet[g2[indexToReturn][0]]+jet[g2[indexToReturn][1]]+lepton[0]+lepton[1];
          LorentzTransform LT=LorentzTransform::mkFrameTransform(CM);
          p_T_ISR = (LT.transform(jet[g1[indexToReturn][0]]+jet[g1[indexToReturn][1]])).pT();
          p_T_I = (LT.transform(EtMiss1)).pT();
          p_T = (CM).pT();

          Vector3 p_I = (LT.transform(EtMiss1)).p3();
          Vector3 p_S = (LT.transform(EtMiss1+jet[g2[indexToReturn][0]]+jet[g2[indexToReturn][1]]+lepton[0]+lepton[1])).p3();
          R_ISR = (p_S.dot(p_I))/(p_S.dot(p_S));
          delta_phi = deltaPhi((LT.transform(EtMiss1)).p3(),(LT.transform(jet[g1[indexToReturn][0]]+jet[g1[indexToReturn][1]])).p3());
        }

        if (delta_phi < 2.8) break;
        _cutflow2lISR.fill(5);

        if (R_ISR < 0.4 || R_ISR > 0.75) break;
        _cutflow2lISR.fill(6);

        if (p_T_ISR < 180*GeV) break;
        _cutflow2lISR.fill(7);

        if (p_T_I < 100*GeV) break;
        _cutflow2lISR.fill(8);

        if (p_T > 20*GeV) break;
        _cutflow2lISR.fill(9);

        break;
      }


      // 3-lepton High (n=0), Int (n=1) and Low (n=2) Selection
      for (int n = 0; n < 3; ++n) {

        while (true) {

          // Require the leading electron and the leading muon to have pT > 25 GeV
          if (elecs.size() != 0 && elecs[0].pT() < 25*GeV) break;
          if (muons.size() != 0 && muons[0].pT() < 25*GeV) break;

          // Require 3 leptons and a pair of leptons with opposite charge and same flavour
          double n_leptons = leptons.size();
          if (n_leptons != 3) break;

          vector<vector<double>> g = {{1,2},{0,2},{0,1}};
          int indexToReturn = 3;
          int indexValue = 100000;
          int newValue = 0;
          for (int i = 0; i < 3; i++) {
            if (!(leptons[g[i][0]].abspid() == leptons[g[i][1]].abspid() && leptons[g[i][0]].charge() != leptons[g[i][1]].charge())) continue;
            newValue = abs((leptons[g[i][0]].mom()+leptons[g[i][1]].mom()).mass() - 90*GeV);
            if (newValue <= indexValue) {
              indexToReturn = i;
              indexValue = newValue;
            }
          }

          if (indexToReturn == 3) break;

          _cutflow3l[n].fill(1);

          // Requirements on n_jets
          double n_jets = jets.size();
          if (n != 2) {
            if (n_jets >= 3) break;
            if (any(jets, hasBTag())) break;
          } else {
            if (n_jets != 0) break;
          }

          // Obtain 4-momentum of the leptons and implement the requirements on their transverse momentum
          vector<FourMomentum> lepton;
          for (int l = 0; l < 3; ++l) {
            lepton.push_back(leptons[l].mom());
          }

          double p_l1_T = lepton[0].pT();
          if (p_l1_T < 60) break;

          if (n==0) {
            double p_l2_T = lepton[1].pT();
            if (p_l2_T < 60*GeV) break;
          } else if (n==1) {
            double p_l2_T = lepton[1].pT();
            if (p_l2_T < 50*GeV) break;
          } else {
            double p_l2_T = lepton[1].pT();
            if (p_l2_T < 40*GeV) break;
          }

          if (n==0) {
            double p_l3_T = (lepton[2]).pT();
            if (p_l3_T < 40*GeV) break;
          } else {
            double p_l3_T = (lepton[2]).pT();
            if (p_l3_T < 30*GeV) break;
          }

          // Pre-selection Cut
          _cutflow3l[n].fill(2);

          // Requirement on m_ll
          double m_ll = (lepton[g[indexToReturn][0]]+lepton[g[indexToReturn][1]]).mass();
          if (m_ll < 75*GeV || m_ll > 105*GeV) break;

          _cutflow3l[n].fill(3);

          // Obtain the missing transverse momentum vector
          Vector3 EtMissX = apply<SmearedMET>(event,"MET").vectorPt();
          Vector3 EtMiss  = -EtMissX;

          // Requirement on m_W_T
          double deltaphi = deltaPhi(EtMiss,(lepton[indexToReturn]).p3());
          double m_W_T = sqrt(2*((lepton[indexToReturn]).pT())*sqrt(EtMiss.dot(EtMiss))*(1-cos(deltaphi)));

          if (n==0) {
            if (m_W_T < 150*GeV) break;
          } else if (n==1) {
            if (m_W_T < 130*GeV) break;
          } else {
            if (m_W_T < 100*GeV) break;
          }

          _cutflow3l[n].fill(4);

          // Invisible mass JR and Invisible Rapidity JR to obtain Invisible System 4-momentum
          FourMomentum P_V = lepton[0]+lepton[1]+lepton[2];
          double M_I = sqrt(P_V.mass2() - 4*m_ll*((lepton[indexToReturn]).mass()));
          double M_V = P_V.mass();
          double p_I_z = (P_V.pz())*sqrt(EtMiss.dot(EtMiss)+sqr(M_I)) / sqrt(P_V.pT2() + sqr(M_V));

          FourMomentum P_I;
          P_I.setPM(EtMiss.x(), EtMiss.y(), p_I_z, M_I);

          // Lorentz boost to CM frame
          LorentzTransform LT = LorentzTransform::mkFrameTransform(P_I+P_V);
          FourMomentum P_F_I; //4-momentum of invisible system in CM frame
          FourMomentum P_F_Va; //4-momentum of visible system a in CM frame
          FourMomentum P_F_Vb; //4-momentum of visible system b in CM frame
          FourMomentum P_F_V; //4-momentum of visible system in CM frame
          FourMomentum P_F_Va1; //4-momentum of visible system a1 in CM frame
          FourMomentum P_F_Va2; //4-momentum of visible system a2 in CM frame
          FourMomentum P_F_Vb1; //4-momentum of visible system b1 in CM frame
          FourMomentum P_F_Vb2; //4-momentum of visible system b2 in CM frame
          if ((lepton[indexToReturn]).mass() > (lepton[g[indexToReturn][0]]+lepton[g[indexToReturn][1]]).mass()) {
            P_F_I = LT.transform(P_I);
            P_F_Va = LT.transform(lepton[indexToReturn]);
            P_F_Vb = LT.transform(lepton[g[indexToReturn][0]]+lepton[g[indexToReturn][1]]);
            P_F_V = LT.transform(P_V);
            //
            P_F_Va1 = LT.transform(lepton[indexToReturn]);
            P_F_Vb1 = LT.transform(lepton[g[indexToReturn][0]]);
            P_F_Vb2 = LT.transform(lepton[g[indexToReturn][1]]);

          } else {

            P_F_I = LT.transform(P_I);
            P_F_Va = LT.transform(lepton[g[indexToReturn][0]]+lepton[g[indexToReturn][1]]);
            P_F_Vb = LT.transform(lepton[indexToReturn]);
            P_F_V = LT.transform(P_V);
            //
            P_F_Va1 = LT.transform(lepton[g[indexToReturn][0]]);
            P_F_Va2 = LT.transform(lepton[g[indexToReturn][1]]);
            P_F_Vb1 = LT.transform(lepton[indexToReturn]);
          }

          // Obtain variables defined in the Contra-boost Invariant JR
          double M2c = 2*((P_F_Va.E())*(P_F_Vb.E())+(P_F_Va.p3()).dot(P_F_Vb.p3()));
          double m2Va = sqr(P_F_Va.mass());
          double m2Vb = sqr(P_F_Vb.mass());
          double ka = m2Va-m2Vb+M2c-2*sqrt(m2Va)*sqrt(m2Vb);
          double kb = m2Vb-m2Va+M2c-2*sqrt(m2Va)*sqrt(m2Vb);
          double kn = ka*m2Va - kb*m2Vb + M2c*(kb-ka)/2 + (0.5)*sqrt(pow(ka+kb,2)*(pow(M2c,2)-4*m2Va*m2Vb));
          double k2d = sqr(ka)*m2Va+sqr(kb)*m2Vb+ka*kb*M2c;
          double k = kn/k2d;
          double ca = (1+k*ka)/2;
          double cb = (1+k*kb)/2;
          double c = 0.5*(P_F_V.E()+sqrt(pow(P_F_V.E(),2)+pow(M_I,2)-pow(P_V.mass(),2)))/(ca*(P_F_Va.E())+cb*(P_F_Vb.E()));

          // Apply Contra-boost Invariant JR to obtain invisible particles' 4-momenta
          double p_F_Iax = (P_F_Va.px())*(c*ca-1)-(P_F_Vb.px())*c*cb;
          double p_F_Iay = (P_F_Va.py())*(c*ca-1)-(P_F_Vb.py())*c*cb;
          double p_F_Iaz = (P_F_Va.pz())*(c*ca-1)-(P_F_Vb.pz())*c*cb;
          double p_F_Ibx = (P_F_Vb.px())*(c*cb-1)-(P_F_Va.px())*c*ca;
          double p_F_Iby = (P_F_Vb.py())*(c*cb-1)-(P_F_Va.py())*c*ca;
          double p_F_Ibz = (P_F_Vb.pz())*(c*cb-1)-(P_F_Va.pz())*c*ca;
          double E_F_Ia = (c*ca-1)*(P_F_Va.E())+c*cb*(P_F_Vb.E());
          double E_F_Ib = (c*cb-1)*(P_F_Vb.E())+c*ca*(P_F_Va.E());

          FourMomentum P_F_Ia;
          FourMomentum P_F_Ib;
          P_F_Ia.setPE(p_F_Iax,p_F_Iay,p_F_Iaz,E_F_Ia);
          P_F_Ib.setPE(p_F_Ibx,p_F_Iby,p_F_Ibz,E_F_Ib);

          // Lorentz Transform to the Pb frame from the CM frame
          LorentzTransform LTP;
          if ((lepton[indexToReturn]).mass() > (lepton[g[indexToReturn][0]]+lepton[g[indexToReturn][1]]).mass()) {
            LTP = LorentzTransform::mkFrameTransform(P_F_Vb+P_F_Ib);
          } else {
            LTP = LorentzTransform::mkFrameTransform(P_F_Va+P_F_Ia);
          }

          // p_lab_T_PP/(p_lab_T_PP+p_lab_T_31) requirement
          LorentzTransform LTR = LorentzTransform::mkObjTransform(P_I+P_V);
          double p_lab_T_PP = (LTR.transform(P_F_Va+P_F_Vb+P_F_Ia+P_F_Ib)).pT();
          double H_PP_T_31;
          double H_PP_31;
          double H_Pb_11;
          double H_Pb_21;
          if ((lepton[indexToReturn]).mass() > (lepton[g[indexToReturn][0]]+lepton[g[indexToReturn][1]]).mass()) {
            H_PP_T_31 = P_F_Va1.pT()+P_F_Vb1.pT()+P_F_Vb2.pT()+(P_F_Ia+P_F_Ib).pT();
            H_PP_31 = P_F_Va1.p()+P_F_Vb1.p()+P_F_Vb2.p()+(P_F_Ia+P_F_Ib).p();
            H_Pb_11 = (LTP.transform(P_F_Vb1+P_F_Vb2)).p()+(LTP.transform(P_F_Ib)).p();
            H_Pb_21 = (LTP.transform(P_F_Vb1)).p()+(LTP.transform(P_F_Vb2)).p()+(LTP.transform(P_F_Ib)).p();
          } else {
            H_PP_T_31 = P_F_Va1.pT()+P_F_Va2.pT()+P_F_Vb1.pT()+(P_F_Ia+P_F_Ib).pT();
            H_PP_31 = P_F_Va1.p()+P_F_Va2.p()+P_F_Vb1.p()+(P_F_Ia+P_F_Ib).p();
            H_Pb_11 = (LTP.transform(P_F_Va1+P_F_Va2)).p()+(LTP.transform(P_F_Ia)).p();
            H_Pb_21 = (LTP.transform(P_F_Va1)).p()+(LTP.transform(P_F_Va2)).p()+(LTP.transform(P_F_Ia)).p();
          }

          double fraction1 = H_PP_T_31/H_PP_31;
          if (n==0) {
            if (fraction1 < 0.75) break;
          } else if (n==1) {
            if (fraction1 < 0.8) break;
          } else {
            if (fraction1 < 0.9) break;
          }
          _cutflow3l[n].fill(5);

          double fraction2 = H_Pb_11/H_Pb_21;
          if (n==0) {
            if (fraction2 < 0.8) break;
            _cutflow3l[n].fill(6);
          } else if (n==1) {
            if (fraction2 < 0.75) break;
            _cutflow3l[n].fill(6);
          } else { /* ??? */ }

          if (n==0) {
            if (H_PP_31 < 550*GeV) break;
            _cutflow3l[n].fill(7);
          } else if (n==1) {
            if (H_PP_31 < 450*GeV) break;
            _cutflow3l[n].fill(7);
          } else {
            if (H_PP_31 < 250*GeV) break;
            _cutflow3l[n].fill(6);
          }

          double fraction3 = p_lab_T_PP/(p_lab_T_PP+H_PP_T_31);
          if (n==0) {
            if (fraction3 > 0.2) break;
            _cutflow3l[n].fill(8);
          } else if (n==1) {
            if (fraction3 > 0.15) break;
            _cutflow3l[n].fill(8);
          } else {
            if (fraction3 > 0.05) break;
            _cutflow3l[n].fill(7);
          }

          break;
        }
      }


      // 3-lepton ISR Selection
      while (true) {

        // Require the leading electron and the leading muon to have pT > 25GeV
        if (elecs.size() != 0 && elecs[0].pT() < 25*GeV) break;
        if (muons.size() != 0 && muons[0].pT() < 25*GeV) break;

        // Require only 3 leptons
        double n_leptons = leptons.size();
        if (n_leptons != 3) break;

        // Require and identify the pair of opposite-charged and same-flavoured leptons
        vector<vector<double>> g = {{1,2},{0,2},{0,1}};
        int indexToReturn = 3;
        int indexValue = 100000;
        int newValue = 0;
        for (int i = 0; i < 3; i++) {
          if (!(leptons[g[i][0]].abspid() == leptons[g[i][1]].abspid() && leptons[g[i][0]].charge() != leptons[g[i][1]].charge())) continue;
          newValue = abs(((leptons[g[i][0]]).mom()+(leptons[g[i][1]]).mom()).mass()-90);
          if (newValue <= indexValue) {
            indexToReturn = i;
            indexValue = newValue;
          }
        }

        if (indexToReturn == 3) break;
        _cutflow3lISR.fill(1);

        // Requirements on jet number and forbid B-tags
        double n_jets=jets.size();
        if (n_jets < 1 || n_jets > 3) break;
        if (any(jets, hasBTag())) break;

        // Obtain the 4-momenta of leptons and implement requirements on them
        vector<FourMomentum> lepton;
        double M;
        for (int n = 0; n < 3; ++n) {
          lepton.push_back(leptons[n].mom());
          M = lepton[n].mass();
          lepton[n].setPz(0);
          lepton[n].setE(sqrt(pow(lepton[n].px(),2)+pow(lepton[n].py(),2)+pow(M,2)));
        }
        if (lepton[0].pT() < 25*GeV || lepton[1].pT() < 25*GeV || lepton[2].pT() < 20*GeV) break;

        // Pre-selection Cut
        _cutflow3lISR.fill(2);

        // Obtain the missing momentum 3-vector
        Vector3 EtMiss = apply<SmearedMET>(event,"MET").vectorMissingPt();

        // Requirement on m_ll
        double m_ll = (leptons[g[indexToReturn][0]].mom()+leptons[g[indexToReturn][1]].mom()).mass();
        if (m_ll < 75*GeV || m_ll > 105*GeV) break;
        _cutflow3lISR.fill(3);

        // Requirement on m_W_T
        double deltaphi = deltaPhi(EtMiss,(lepton[indexToReturn]).p3());
        double m_W_T = sqrt(2*((lepton[indexToReturn]).pT())*sqrt(EtMiss.dot(EtMiss))*(1-cos(deltaphi)));
        if (m_W_T < 100*GeV) break;
        _cutflow3lISR.fill(4);

        // Obtain the missing transverse 4-momentum
        FourMomentum EtMiss1;
        EtMiss1.setPM(EtMiss.x(),EtMiss.y(),EtMiss.z(),sqrt(EtMiss.dot(EtMiss)));
        double p_T_ISR;
        double p_T_I;
        double p_T;
        double R_ISR;
        double delta_phi;
        FourMomentum ISR;

        // Obtain the 4-momentum of the ISR system
        vector<FourMomentum> jet;
        for (int n = 0; n < n_jets; ++n) {
          jet.push_back(jets[n].momentum());
          M = jet[n].mass();
          jet[n].setPz(0);
          jet[n].setE(sqrt(pow(jet[n].px(),2)+pow(jet[n].py(),2)+pow(M,2)));
        }

        if (n_jets==1) {
          ISR = jet[0];
        } else if (n_jets==2) {
          ISR = jet[0]+jet[1];
        } else {
          ISR = jet[0]+jet[1]+jet[2];
        }

        // Four Momentum of the CM System
        FourMomentum CM = EtMiss1+ISR+lepton[0]+lepton[1]+lepton[2];
        LorentzTransform LT = LorentzTransform::mkFrameTransform(CM);
        p_T_ISR = (LT.transform(ISR)).pT();
        p_T_I = (LT.transform(EtMiss1)).pT();
        p_T = (CM).pT();

        Vector3 p_I = (LT.transform(EtMiss1)).p3();
        Vector3 p_S = (LT.transform(EtMiss1+lepton[0]+lepton[1]+lepton[2])).p3();
        R_ISR = (p_S.dot(p_I))/(p_S.dot(p_S));
        delta_phi = angle((LT.transform(EtMiss1)).p3(),(LT.transform(ISR)).p3());

        if (delta_phi < 2.0) break;
        _cutflow3lISR.fill(5);

        if (R_ISR < 0.55 || R_ISR > 1.0) break;
        _cutflow3lISR.fill(6);

        if (p_T_ISR < 100*GeV) break;
        _cutflow3lISR.fill(7);

        if (p_T_I < 80*GeV) break;
        _cutflow3lISR.fill(8);

        if (p_T > 25*GeV) break;
        _cutflow3lISR.fill(9);

        break;
      }

    }


    /// Finalise cutflow scaling etc.
    void finalize() {

      _cutflow2l[0].normalize(1673, 0);
      _cutflow2l[1].normalize(4369, 0);
      _cutflow2l[2].normalize(65247, 0);
      _cutflow2lISR.normalize(65247, 0);
      _cutflow3l[0].normalize(1673, 0);
      _cutflow3l[1].normalize(4369, 0);
      _cutflow3l[2].normalize(65247, 0);
      _cutflow3lISR.normalize(65247, 0);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow2l[0]);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow2l[1]);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow2l[2]);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow2lISR);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow3l[0]);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow3l[1]);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow3l[2]);
      MSG_INFO("CUTFLOWS:\n\n" << _cutflow3lISR);
    }

    Cutflows _cutflow2lHigh;
    Cutflows _cutflow2lInt;
    Cutflows _cutflow2lLow;
    Cutflows _cutflow2lISR;
    Cutflows _cutflow3lHigh;
    Cutflows _cutflow3lInt;
    Cutflows _cutflow3lLow;
    Cutflows _cutflow3lISR;

    vector<Cutflows> _cutflow2l={_cutflow2lHigh,_cutflow2lInt,_cutflow2lLow};
    vector<Cutflows> _cutflow3l={_cutflow3lHigh,_cutflow3lInt,_cutflow3lLow};

  };


  RIVET_DECLARE_PLUGIN(ATLAS_2018_I1676551);

}
