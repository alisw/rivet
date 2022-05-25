// -*- C++ -*-

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"

#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"


namespace Rivet {


  /// @brief Jet substructure at 13 TeV
  class ATLAS_2019_I1724098: public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2019_I1724098);

  private:

    /// @name Analysis methods
    //@{

    void init() {

      _mode = 0; // default is do eveything
      if ( getOption("MODE") == "DJ" )       _mode = 1;
      else if ( getOption("MODE") == "TW" )  _mode = 2;

      //Projections
      const FinalState fs(Cuts::abseta < 4.5);

      // Get photons used to dress leptons
      const FinalState photons(Cuts::abspid == PID::PHOTON);

      // Use all bare muons as input to the DressedMuons projection
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true);
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true);

      // Muons must have |eta| < 2.5
      Cut eta_ranges = Cuts::abseta < 2.5;
      DressedLeptons dressed_mu(photons, bare_mu, 0.1, eta_ranges && Cuts::pT > 30*GeV, true);
      declare(dressed_mu, "muons");
      DressedLeptons dressed_el(photons, bare_el, 0.1, eta_ranges && Cuts::pT > 25*GeV, true);
      declare(dressed_el, "electrons");

      FastJets fj(fs, FastJets::ANTIKT, 1.0, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(fj, "FJets");
      FastJets sj(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      declare(sj, "Jets");

      ChargedLeptons lfs(FinalState(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV));
      declare(lfs, "LFS");

      MissingMomentum missmom(fs);
      declare(missmom, "MissingMomentum");

      _trimmer = fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05));

      // Dijet
      if (_mode == 0 || _mode == 1) {
        // SD
        book(_h["dj_sdnsj"],1,1,1);
        book(_h["dj_sdlha"]  ,2,1,1);
        book(_h["dj_sdc2"]   ,3,1,1);
        book(_h["dj_sdd2"]   ,4,1,1);
        book(_h["dj_sdecf2"] ,5,1,1);
        book(_h["dj_sdecf3"] ,6,1,1);
        // TRIMMED
        book(_h["dj_nsj"]  ,23,1,1);
        book(_h["dj_lha"]  ,24,1,1);
        book(_h["dj_c2"]   ,25,1,1);
        book(_h["dj_d2"]   ,26,1,1);
        book(_h["dj_ecf2"] ,27,1,1);
        book(_h["dj_ecf3"] ,28,1,1);
      }

      if (_mode == 0 || _mode == 2) {
        //Top SD
        book(_h["tw_sdnsj"]   ,7,1,1);
        book(_h["tw_sdlha"]   ,8,1,1);
        book(_h["tw_sdc2"]    ,9,1,1);
        book(_h["tw_sdd2"]    ,10,1,1);
        book(_h["tw_sdecf2"]  ,11,1,1);
        book(_h["tw_sdecf3"]  ,12,1,1);
        book(_h["tw_sdtau21"] ,13,1,1);
        book(_h["tw_sdtau32"] ,14,1,1);
         //W SD
        book(_h["tw_wsdnsj"]   ,15,1,1);
        book(_h["tw_wsdlha"]   ,16,1,1);
        book(_h["tw_wsdc2"]    ,17,1,1);
        book(_h["tw_wsdd2"]    ,18,1,1);
        book(_h["tw_wsdecf2"]  ,19,1,1);
        book(_h["tw_wsdecf3"]  ,20,1,1);
        book(_h["tw_wsdtau21"] ,21,1,1);
        book(_h["tw_wsdtau32"] ,22,1,1);
        //Top TRIMMED
        book(_h["tw_nsj"]   ,29,1,1);
        book(_h["tw_lha"]   ,30,1,1);
        book(_h["tw_c2"]    ,31,1,1);
        book(_h["tw_d2"]    ,32,1,1);
        book(_h["tw_ecf2"]  ,33,1,1);
        book(_h["tw_ecf3"]  ,34,1,1);
        book(_h["tw_tau21"] ,35,1,1);
        book(_h["tw_tau32"] ,36,1,1);
        //W TRIMMED
        book(_h["tw_wnsj"]   ,37,1,1);
        book(_h["tw_wlha"]   ,38,1,1);
        book(_h["tw_wc2"]    ,39,1,1);
        book(_h["tw_wd2"]    ,40,1,1);
        book(_h["tw_wecf2"]  ,41,1,1);
        book(_h["tw_wecf3"]  ,42,1,1);
        book(_h["tw_wtau21"] ,43,1,1);
        book(_h["tw_wtau32"] ,44,1,1);
      }

    }

    /// Do the analysis
    void analyze(const Event& event) {

      if (_mode == 0 || _mode == 1)  doDIJET(event);
      if (_mode == 0 || _mode == 2)  doTW(event);

    }


    void doDIJET(const Event& event) {

      double  nsub, lha, ecf2, ecf3, c2, d2;
      double  sdnsub, sdlha, sdecf2, sdecf3, sdc2, sdd2;

      lha = sdlha = 0.0;

      const double beta = 1;

      const Particles & leptons = apply<ChargedLeptons>(event, "LFS").particles();
      if (leptons.size())  return;

      // Normal fatjets
      const Jets &fjets = apply<JetAlg>(event, "FJets").jetsByPt();

      // Trim the fatjets
      PseudoJets tr_ljets;
      for (size_t i = 0; i < fjets.size(); ++i) {
        tr_ljets += _trimmer(fjets[i]);
        tr_ljets[tr_ljets.size()-1].set_user_index(i);
      }

      size_t nBaseline = count(tr_ljets, [](const Jet &j) { return j.pT() > 200*GeV && j.abseta() < 2.5; });
      if (nBaseline < 2)  return;

      ifilter_select(tr_ljets, [](const PseudoJet &j) { return j.perp() > 450*GeV; });
      if (tr_ljets.size() > 1)  tr_ljets = sorted_by_pt(tr_ljets);
      else if (tr_ljets.empty())  return;

      if (abs(tr_ljets[0].eta()) > 1.5)  return;
      const fastjet::PseudoJet &LJet = tr_ljets[0];
      size_t uindex = tr_ljets[0].user_index();

      // Nsubjets
      JetDefinition subjet_def(fastjet::kt_algorithm, 0.2);
      ClusterSequence subjet_cs(LJet.constituents(), subjet_def);
      PseudoJets subjets = sorted_by_pt(subjet_cs.inclusive_jets(10.0));
      nsub = subjets.size();

      // LHA
      for (const PseudoJet& p : LJet.constituents()) {
        double trpt = p.pt();
        double trtheta = p.squared_distance(LJet);
        lha += pow(trpt, 1.0) * pow(trtheta, 0.25);
      }
      double lterm = pow(LJet.pt(), 1.0) * pow(1.0, 0.5);
      if (lterm !=0)  lha /= lterm;
      else            lha = -99;

      // C2
      fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
      fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
      fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);
      double recf3 = ECF3(LJet);
      double recf2 = ECF2(LJet);
      double recf1 = ECF1(LJet);
      c2 = (recf2 != 0 ? recf3 * recf1 / (recf2*recf2) : -1);
      d2 = (recf2 != 0 ? recf3 * (recf1*recf1*recf1) /(recf2*recf2*recf2) : -1);
      ecf2 = (recf1 != 0 ? recf2 / (recf1*recf1) : -1);
      ecf3 = (recf1 != 0 ? recf3 / (recf1*recf1*recf1) : -1);

      // Fill Histograms for trimmed
      _h["dj_nsj"]->fill(nsub);
      _h["dj_c2"]->fill(c2);
      _h["dj_d2"]->fill(d2);
      _h["dj_lha"]->fill(lha);
      _h["dj_ecf2"]->fill(ecf2);
      _h["dj_ecf3"]->fill(ecf3);


      ////////////////////////////////////////////
      // Soft Drop
      fastjet::contrib::SoftDrop sd(0.0, 0.1);
      PseudoJet SDLJet = sd(fjets[uindex]);

      ClusterSequence subjet_sdcs(SDLJet.constituents(), subjet_def);

      PseudoJets sdsubjets = sorted_by_pt(subjet_sdcs.inclusive_jets(10.0));
      sdnsub = sdsubjets.size();

      for (const PseudoJet& sd_p : SDLJet.constituents()) {
        double spt = sd_p.pt();
        double stheta = sd_p.squared_distance(SDLJet);
        sdlha += pow(spt, 1.0) * pow(stheta, 0.25);
      }

      double sdterm = pow(SDLJet.pt(), 1.0) * pow(1.0, 0.5);
      if (sdterm !=0)  sdlha /= sdterm;
      else             sdlha = -99;

      double sdrecf3 = ECF3(SDLJet);
      double sdrecf2 = ECF2(SDLJet);
      double sdrecf1 = ECF1(SDLJet);

      sdc2 = (sdrecf2 != 0 ? sdrecf3 * sdrecf1 / (sdrecf2*sdrecf2) : -1);
      sdd2 = (sdrecf2 != 0 ? sdrecf3 * (sdrecf1*sdrecf1*sdrecf1) /(sdrecf2*sdrecf2*sdrecf2) : -1);

      sdecf2 = (sdrecf1 !=0 ? sdrecf2 /(sdrecf1*sdrecf1) : -1);
      sdecf3 = (sdrecf1 !=0 ? sdrecf3 / (sdrecf1*sdrecf1*sdrecf1) : -1);

      _h["dj_sdnsj"]->fill(sdnsub);
      _h["dj_sdc2"]->fill(sdc2);
      _h["dj_sdd2"]->fill(sdd2);
      _h["dj_sdlha"]->fill(sdlha);
      _h["dj_sdecf2"]->fill(sdecf2);
      _h["dj_sdecf3"]->fill(sdecf3);

    }

    void doTW(const Event& event) {

      double nsub, lha, tau21, tau32, ecf2, ecf3, c2, d2;
      double sdnsub, sdlha, sdtau21, sdtau32, sdecf2, sdecf3, sdc2, sdd2;
      double wnsub, wlha, wtau21, wtau32, wecf2, wecf3, wc2, wd2;
      double wsdnsub, wsdlha, wsdtau21, wsdtau32, wsdecf2, wsdecf3, wsdc2, wsdd2;

      lha = sdlha = wlha = wsdlha = 0.0;

      const double beta = 1, Rcut = 1;

      const vector<DressedLepton>& muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
      if (muons.size() != 1)  return;

      const vector<DressedLepton>& electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();
      if (electrons.size() != 0)  return;

      const FourMomentum muonmom = muons[0].momentum();
      const MissingMomentum& missmom = apply<MissingMomentum>(event, "MissingMomentum");
      FourMomentum missvec = missmom.missingMomentum();
      double met =  missmom.missingPt();
      if (met < 20*GeV)  return;
      const double transmass = sqrt( 2 * muons[0].pT() * met * (1 - cos(deltaPhi(muons[0], missvec))) );
      if (transmass + met <= 60*GeV)  return;

      const Jets& jets  = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      if(jets.empty())  return;

      int lepJetIndex = -1;
      for (size_t i = 0; i < jets.size(); ++i) {
        const Jet& jet = jets[i];
        if ((deltaR(jet, muons[0]) < 1.5) && (deltaR(jet, muons[0]) > 0.4) ) {
          lepJetIndex = i;
          break;
        }
      }
      if (lepJetIndex < 0)  return;
      const Jet& lepjet = jets[lepJetIndex];


      const Jets fjets = apply<JetAlg>(event, "FJets").jetsByPt();
      PseudoJets tr_ljets_all;
      for (const Jet& j : fjets) {
        tr_ljets_all += _trimmer(j);
      }

      PseudoJets tr_ljets;
      for (size_t i = 0; i < tr_ljets_all.size(); ++i) {
        const PseudoJet tj = tr_ljets_all[i];
        if (tj.perp() > 150*GeV && fabs(tj.eta()) < 2.5) {
          tr_ljets += tj;
          tr_ljets[tr_ljets.size()-1].set_user_index(i);
        }
      }

      if (tr_ljets.size() < 1)  return;
      if (tr_ljets.size() > 1)  tr_ljets = sorted_by_pt(tr_ljets);

      const Jet& fjet = tr_ljets[0];
      size_t uindex = tr_ljets[0].user_index();

      const double dR_fatjet = deltaR(lepjet, fjet);
      const double dPhi_fatjet = deltaPhi(muons[0], fjet);

      if (dR_fatjet < 1.5 || dPhi_fatjet < 2.3)  return;

      Jets bjets, non_bjets;
      for (const Jet& jet : jets)
        (jet.bTagged() ? bjets : non_bjets) += jet;
      if (bjets.empty())  return;

      double min_bdR = 99;
      int bindex = 0;
      int k=0;

      for (const Jet& bjet : bjets) {
        double bdR = deltaR(fjet, bjet);
        if(bdR <  min_bdR){
	        min_bdR = bdR;
          bindex = k;
        }
    	  k++;
      }

      size_t tw = 0;
      double dR_fjet_bjet = deltaR(fjet, bjets[bindex]);
      if (fjet.abseta() < 1.5) {
        if (dR_fjet_bjet < 1.0 && fjet.mass() > 140 && fjet.pT() > 350*GeV)  tw = 1;
        if (dR_fjet_bjet > 1.0 && dR_fjet_bjet < 1.8 && fjet.mass() < 100*GeV && fjet.mass() > 60*GeV && fjet.pT() > 200*GeV) tw = 2;
      }


      // Top plots:
      if (tw==1) {

        PseudoJet LJet = fjet;
        JetDefinition subjet_def(fastjet::kt_algorithm, 0.2);
        ClusterSequence subjet_cs(LJet.constituents(), subjet_def);
        PseudoJets subjets = sorted_by_pt(subjet_cs.inclusive_jets(10.0));
        nsub = subjets.size();

        // LHA
        for (const PseudoJet& p : LJet.constituents()){
          double pt = p.pt();
          double theta = p.squared_distance(LJet);
          lha += pow(pt, 1.0) * pow(theta, 0.25);
        }
        double lterm = pow(LJet.pt(), 1.0) * pow(1.0, 0.5);
        if (lterm)  lha /= lterm;
        else        lha = -99;


        // NSubjettiness
        fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        double tau1 = nSub1.result(LJet);
        double tau2 = nSub2.result(LJet);
        double tau3 = nSub3.result(LJet);
        if(tau1 != 0) tau21 = tau2/tau1;
        else tau21 = -99;
        if(tau2 != 0) tau32 = tau3/tau2;
        else tau32 = -99;


        //C2
        fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);

        double recf3 = ECF3(LJet);
        double recf2 = ECF2(LJet);
        double recf1 = ECF1(LJet);


        c2 = (recf2 != 0 ? recf3 * recf1 / (recf2*recf2) : -1);
        d2 = (recf2 != 0 ? recf3 * (recf1*recf1*recf1) /(recf2*recf2*recf2) : -1);

        ecf2 = (recf1 !=0 ? recf2 /(recf1*recf1) : -1);
        ecf3 = (recf1 !=0 ? recf3 / (recf1*recf1*recf1) : -1);

        _h["tw_nsj"]->fill(nsub);
        _h["tw_lha"]->fill(lha);
        _h["tw_tau21"]->fill(tau21);
        _h["tw_tau32"]->fill(tau32);
        _h["tw_c2"]->fill(c2);
        _h["tw_d2"]->fill(d2);
        _h["tw_ecf2"]->fill(ecf2);
        _h["tw_ecf3"]->fill(ecf3);


        // Soft Drop
        fastjet::contrib::SoftDrop sd(0.0, 0.1);
        PseudoJet SDLJet = sd(fjets[uindex]);
        ClusterSequence subjet_sdcs(SDLJet.constituents(), subjet_def);

        PseudoJets sdsubjets = sorted_by_pt(subjet_sdcs.inclusive_jets(10.0));
        sdnsub = sdsubjets.size();

        for (const PseudoJet& sd_p : SDLJet.constituents()){
          double spt = sd_p.pt();
          double stheta = sd_p.squared_distance(SDLJet);
          sdlha += pow(spt, 1.0) * pow(stheta, 0.25);
        }

        double sdlterm = pow(SDLJet.pt(), 1.0) * pow(1.0, 0.5);
        if (sdlterm)  sdlha /= sdlterm;
        else          sdlha = -99;


        double sdtau1 = nSub1.result(SDLJet);
        double sdtau2 = nSub2.result(SDLJet);
        double sdtau3 = nSub3.result(SDLJet);
        if(sdtau1 != 0) sdtau21 = sdtau2/sdtau1;
        else sdtau21 = -99;
        if(sdtau2 != 0) sdtau32 = sdtau3/sdtau2;
        else sdtau32 = -99;

        double sdrecf3 = ECF3(SDLJet);
        double sdrecf2 = ECF2(SDLJet);
        double sdrecf1 = ECF1(SDLJet);

        sdc2 = (sdrecf2 != 0 ? sdrecf3 * sdrecf1 / (sdrecf2*sdrecf2) : -1);
        sdd2 = (sdrecf2 != 0 ? sdrecf3 * (sdrecf1*sdrecf1*sdrecf1) /(sdrecf2*sdrecf2*sdrecf2) : -1);

        sdecf2 = (sdrecf1 !=0 ? sdrecf2 /(sdrecf1*sdrecf1) : -1);
        sdecf3 = (sdrecf1 !=0 ? sdrecf3 / (sdrecf1*sdrecf1*sdrecf1) : -1);

        _h["tw_sdnsj"]->fill(sdnsub);
        _h["tw_sdlha"]->fill(sdlha);
        _h["tw_sdtau21"]->fill(sdtau21);
        _h["tw_sdtau32"]->fill(sdtau32);
        _h["tw_sdc2"]->fill(sdc2);
        _h["tw_sdd2"]->fill(sdd2);
        _h["tw_sdecf2"]->fill(sdecf2);
        _h["tw_sdecf3"]->fill(sdecf3);

      }


      // W plots

      if(tw ==2){

        PseudoJet LJet = fjet;
        JetDefinition subjet_def(fastjet::kt_algorithm, 0.2);
        ClusterSequence subjet_cs(LJet.constituents(), subjet_def);

        PseudoJets subjets = sorted_by_pt(subjet_cs.inclusive_jets(10.0));
        wnsub = subjets.size();

        // LHA
        for (const PseudoJet& wp : LJet.constituents()){
          double wpt = wp.pt();
          double wtheta = wp.squared_distance(fjet);
          wlha += pow(wpt, 1.0) * pow(wtheta, 0.25);
        }

        double fterm = pow(fjet.pt(), 1.0) * pow(1.0, 0.5);
        if (fterm)  wlha /= fterm;
        else        wlha = -99;



        // NSubjettiness

        fastjet::contrib::Nsubjettiness nSub1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        fastjet::contrib::Nsubjettiness nSub2(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        fastjet::contrib::Nsubjettiness nSub3(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::NormalizedMeasure(beta,Rcut));
        double wtau1 = nSub1.result(LJet);
        double wtau2 = nSub2.result(LJet);
        double wtau3 = nSub3.result(LJet);
        if(wtau1 != 0) wtau21 = wtau2/wtau1;
        else wtau21 = -99;
        if(wtau2 != 0) wtau32 = wtau3/wtau2;
        else wtau32 = -99;

        //C2
        fastjet::contrib::EnergyCorrelator ECF3(3,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF2(2,beta,fastjet::contrib::EnergyCorrelator::pt_R);
        fastjet::contrib::EnergyCorrelator ECF1(1,beta,fastjet::contrib::EnergyCorrelator::pt_R);

        double wrecf3 = ECF3(LJet);
        double wrecf2 = ECF2(LJet);
        double wrecf1 = ECF1(LJet);

        wc2 = (wrecf2 != 0 ? wrecf3 * wrecf1 / (wrecf2*wrecf2) : -1);
        wd2 = (wrecf2 != 0 ? wrecf3 * (wrecf1*wrecf1*wrecf1) /(wrecf2*wrecf2*wrecf2) : -1);

        wecf2 = (wrecf1 !=0 ? wrecf2 /(wrecf1*wrecf1) : -1);
        wecf3 = (wrecf1 !=0 ? wrecf3 / (wrecf1*wrecf1*wrecf1) : -1);

        _h["tw_wnsj"]->fill(wnsub);
        _h["tw_wlha"]->fill(wlha);
        _h["tw_wtau21"]->fill(wtau21);
        _h["tw_wtau32"]->fill(wtau32);
        _h["tw_wc2"]->fill(wc2);
        _h["tw_wd2"]->fill(wd2);
        _h["tw_wecf2"]->fill(wecf2);
        _h["tw_wecf3"]->fill(wecf3);


        //SD
        fastjet::contrib::SoftDrop sd(0.0, 0.1);
        PseudoJet SDLJet = sd(fjets[uindex]);
        ClusterSequence subjet_sdcs(SDLJet.constituents(), subjet_def);

        PseudoJets sdsubjets = sorted_by_pt(subjet_sdcs.inclusive_jets(10.0));
        wsdnsub = sdsubjets.size();

        for (const PseudoJet& sd_p : SDLJet.constituents()){
          double spt = sd_p.pt();
          double stheta = sd_p.squared_distance(SDLJet);
          wsdlha += pow(spt, 1.0) * pow(stheta, 0.25);
        }

        double wsdlterm = pow(SDLJet.pt(), 1.0) * pow(1.0, 0.5);
        if (wsdlterm)  wsdlha /= wsdlterm;
        else          wsdlha = -99;


        double wsdtau1 = nSub1.result(SDLJet);
        double wsdtau2 = nSub2.result(SDLJet);
        double wsdtau3 = nSub3.result(SDLJet);
        if (wsdtau1 != 0) wsdtau21 = wsdtau2/wsdtau1;
        else wsdtau21 = -99;
        if (wsdtau2 != 0) wsdtau32 = wsdtau3/wsdtau2;
        else wsdtau32 = -99;

        double wsdrecf3 = ECF3(SDLJet);
        double wsdrecf2 = ECF2(SDLJet);
        double wsdrecf1 = ECF1(SDLJet);

        wsdc2 = (wsdrecf2 != 0 ? wsdrecf3 * wsdrecf1 / (wsdrecf2*wsdrecf2) : -1);
        wsdd2 = (wsdrecf2 != 0 ? wsdrecf3 * (wsdrecf1*wsdrecf1*wsdrecf1) /(wsdrecf2*wsdrecf2*wsdrecf2) : -1);

        wsdecf2 = (wsdrecf1 !=0 ? wsdrecf2 /(wsdrecf1*wsdrecf1) : -1);
        wsdecf3 = (wsdrecf1 !=0 ? wsdrecf3 / (wsdrecf1*wsdrecf1*wsdrecf1) : -1);

        _h["tw_wsdnsj"]->fill(wsdnsub);
        _h["tw_wsdlha"]->fill(wsdlha);
        _h["tw_wsdtau21"]->fill(wsdtau21);
        _h["tw_wsdtau32"]->fill(wsdtau32);
        _h["tw_wsdc2"]->fill(wsdc2);
        _h["tw_wsdd2"]->fill(wsdd2);
        _h["tw_wsdecf2"]->fill(wsdecf2);
        _h["tw_wsdecf3"]->fill(wsdecf3);

      }
    }


    void finalize() {
      for (auto hist : _h) { normalize(hist.second); }
    }

  private:

    fastjet::Filter _trimmer;
    map<string, Histo1DPtr> _h;

  protected:
    size_t _mode;
  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2019_I1724098);
}
