// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {

  /// @brief H->yy differentual cross-sections at 13 TeV
  class ATLAS_2022_I2023464 : public Analysis {
  public:

    /// Default constructor
    RIVET_DEFAULT_ANALYSIS_CTOR( ATLAS_2022_I2023464 );

    /// @name Analysis methods
    //@{
    void init() {
      // Project all final state particles (for missing ET)
      const FinalState fs ( Cuts::abseta < 5.5 ) ;
      declare(fs , "FS");

      // Project final state particles with pT>1 GeV (for track-based isolation)
      declare(ChargedFinalState(Cuts::pT > 1*GeV ), "CFS") ;

      // Project photons which pass kinematic cuts
      // (pT>25 GeV, |eta|<2.37  excluding crack)
      PromptFinalState kin_fs_ph (Cuts::abspid == PID::PHOTON &&
                                  Cuts::pT > 25*GeV &&
                                  Cuts::abseta < 2.37 &&
                                  ( Cuts::abseta < 1.37 || Cuts::abseta > 1.52 ),
                                  true, true);
      declare(kin_fs_ph, "PFS");

      // All photons (to dress leptons)
      FinalState photons( Cuts::abspid == PID::PHOTON );

      // Create dressed mu projection
      // true in last arg = use also muons from prompt tau decays
      PromptFinalState bare_mu(Cuts::abspid == PID::MUON, true);
      DressedLeptons dressed_mu( photons, bare_mu, 0.1, Cuts::abseta < 2.7, true);
      declare(dressed_mu, "MFS");

      // Create dressed e projection
      PromptFinalState bare_el(Cuts::abspid == PID::ELECTRON, true);
      // true in last arg = use also electrons from prompt tau decays
      Cut fid_el = Cuts::abseta < 2.47 && (Cuts::abseta < 1.37 || Cuts::abseta > 1.52);
      DressedLeptons dressed_el( photons, bare_el, 0.1, fid_el, true);
      declare(dressed_el, "EFS");

      // Create AntiKt4TruthWZJets projection
      VetoedFinalState vfs( FinalState(Cuts::abseta < 4.5) );
      vfs.addVetoOnThisFinalState( dressed_el );
      vfs.addVetoOnThisFinalState( dressed_mu );
      FastJets jets( vfs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jets , "jets");

      declare(InvisibleFinalState(true), "MET");

      // Declare histograms

      // plots in main body
      book(_h["fid_regions"], 1, 1, 1); // inclusive regions
      book(_h["pT_yy"], 2, 1, 1); // photon obs, in baseline region

      book(_h["N_j_30"], 3, 1, 1); // N(jet) and N(b-jet) categories
      book(_h["catXS_nbjet"], 4, 1, 1);
      book(_h["pT_j1_30"], 5, 1, 1); // Observables in 1-jet events
      book(_h["pT_yy_JV_30"],  6, 1, 1); // Jet-veto observables
      book(_h["m_jj_30"],  7, 1, 1); // Observables in 2-jet events
      book(_h["Dphi_j_j_30_signed"],  8, 1, 1);
      book(_h["pT_yy_vs_yAbs_yy"],  9, 1, 1); // 2D xsections
      book(_h["VBF_Dphi_j_j_30_signed"], 10, 1, 1); // Observables in VBF fiducial region

      // plots in appendix
      book(_h["rel_pT_y1"], 21, 1, 1); // photon obs"], in baseline region
      book(_h["rel_pT_y2"], 23, 1, 1);
      book(_h["yAbs_yy"], 25, 1, 1);
      book(_h["m_yyj_30"], 27, 1, 1); // Observables in 1-jet events
      book(_h["pT_yyj_30"], 29, 1, 1);
      book(_h["HT_30"], 31, 1, 1);
      book(_h["maxTau_yyj_30"], 33, 1, 1);
      book(_h["sumTau_yyj_30"], 35, 1, 1);
      book(_h["pT_yy_JV_40"], 37, 1, 1); // Jet-veto observables
      book(_h["pT_yy_JV_50"], 39, 1, 1);
      book(_h["pT_yy_JV_60"], 41, 1, 1);
      book(_h["Dphi_yy_jj_30"], 43, 1, 1); // Observables in 2-jet events
      book(_h["pT_yyjj_30"], 45, 1, 1);
      book(_h["pT_yy_vs_pT_yyj"], 47, 1, 1); // 2D xsections
      book(_h["pT_yy_vs_maxTau_yyj"], 49, 1, 1);
      book(_h["rel_DpT_y_y_vs_rel_sumpT_y_y"], 51, 1, 1);
      book(_h["VBF_abs_Zepp"], 53, 1, 1); // Observables in VBF fiducial region
      book(_h["VBF_pT_yyjj_30"], 55, 1, 1);
      book(_h["VBF_pT_j1_30"], 57, 1, 1);
      book(_h["VBF_pT_j1_30_vs_Dphi_j_j_30_signed"], 59, 1, 1);
    }


    void analyze(const Event& event) {

      // Get charged particles
      const Particles& tracks = apply<ChargedFinalState>(event,"CFS").particles();

      // Save prompt photons passing particle-level isolation requirement
      Particles photons;
      for (const Particle& ph : apply<FinalState>(event, "PFS").particlesByPt()) {
        FourMomentum ET_iso (0, 0, 0, 0);
        for (const auto& track : tracks) {
          if ( deltaR(track, ph) > 0.2)  continue;
          ET_iso += track.mom();
        }
        const bool passIso = ET_iso.Et() / ph.pt() < 0.05;
        if ( !passIso )  continue;
        photons += ph;
      }

      // Diphoton selection: two fiducial photons with
      // 105<m_yy<160 GeV
      // pT/m_yy > 0.35, 0.25
      if (photons.size() < 2)  vetoEvent;

      const FourMomentum y1 = photons[0].mom();
      const FourMomentum y2 = photons[1].mom();
      const FourMomentum yy = y1 + y2;
      const double myy = yy.mass();
      if ( y1.pT() < 0.35*myy || y2.pT() < 0.25*myy )  vetoEvent;
      if ( !inRange(myy/GeV, 105., 160.) )  vetoEvent;

      // Retain prompt electrons with pseudorapidity in acceptance
      // and pT>10 GeV (electrons10) or >15 GeV (electrons15)
      // Do not retain electrons overlapping with photons
      Particles electrons10 = apply<FinalState>(event, "EFS").particlesByPt(Cuts::pT > 10*GeV);
      Particles electrons15 = apply<FinalState>(event, "EFS").particlesByPt(Cuts::pT > 15*GeV);
      idiscardIfAnyDeltaRLess(electrons10, photons, 0.4);
      idiscardIfAnyDeltaRLess(electrons15, photons, 0.4);


      // Retain prompt muons with pseudorapidity in acceptance
      // and pT>10 GeV (muons10) or >15 GeV (muons15)
      Particles muons10 = apply<FinalState>(event, "MFS").particlesByPt(Cuts::pT > 10*GeV);
      Particles muons15 = apply<FinalState>(event, "MFS").particlesByPt(Cuts::pT > 15*GeV);
      idiscardIfAnyDeltaRLess(muons10, photons, 0.4);
      idiscardIfAnyDeltaRLess(muons15, photons, 0.4);

      const size_t nlep10 = electrons10.size() + muons10.size();
      const size_t nlep15 = electrons15.size() + muons15.size();

      // Remove jets overlapping with the fiducial photons and electrons
      // Fill vectors of remaining jets (jets30), and subsets of these jets
      // that are central (|eta|<2.5) and b-jets (bjets)
      Jets jets30, jets30central, bjets;
      for (const Jet& jet : apply<FastJets>(event, "jets").jetsByPt(Cuts::pT>30*GeV && Cuts::absrap<4.4)) {
        if ( any(photons,     DeltaRLess(jet, 0.4)) )  continue;
        if ( any(electrons15, DeltaRLess(jet, 0.2)) )  continue;
        jets30 += jet;
        if (jet.abseta() < 2.5) {
          jets30central += jet;
          if (jet.bTagged())  bjets += jet;
        }
      }
      const size_t njets30 = jets30.size();
      const size_t njets30central = jets30central.size();
      const size_t nbjets  = bjets.size();

      // Calculate missing ET
      Vector3 met_vec;
      for (const Particle& p : apply<InvisibleFinalState>(event, "MET").particles()) {
        met_vec += p.mom().perpVec();
      }
      const double met =  met_vec.mod();

      // Fill the histograms

      // Diphoton Xsections
      const double pT_yy = yy.pt();
      _h["pT_yy"]->fill(pT_yy/GeV);

      const double yAbs_yy = yy.absrap();
      _h["yAbs_yy"]->fill(yAbs_yy);

      const double pTy1 = y1.pt();
      _h["rel_pT_y1"]->fill(pTy1/myy);

      const double pTy2 = y2.pt();
      _h["rel_pT_y2"]->fill(pTy2/myy);

      // pT_yy in bins of |y_yy|
      {
        int binYabs(-1), binpT(-1);
        double binWidth=0.0;
        if (pT_yy/GeV < 45.) {
          binpT = 0;
          binWidth = 45.;
        }
        else if (pT_yy/GeV < 120.) {
          binpT = 1;
          binWidth = 120.-45.;
        }
        else if (pT_yy/GeV < 350.) {
          binpT=2;
          binWidth = 350.-120.;
        }

        if (yAbs_yy<0.5) binYabs=0;
        else if (yAbs_yy<1.0) binYabs=1;
        else if (yAbs_yy<1.5) binYabs=2;
        else if (yAbs_yy<2.5) binYabs=3;

        if (binYabs>=0 && binpT>=0) {
          int bin = binYabs*3 + binpT;
          fillHist("pT_yy_vs_yAbs_yy", bin+1, binWidth );
        }
      }

      // [pT(y1)-pT(y2)]/myy in bins of [pT(y1)+pT(y2)]/myy
      const double pTy1py2 = pTy1 + pTy2;
      const double pTy1my2 = pTy1 - pTy2;
      {
        int bin = -1;
        double binWidth = 0.0;
        if ( pTy1py2/myy >= 0.6 && pTy1py2/myy < 0.8) {
          if ( pTy1my2/myy < 0.3 ) {
            bin = 0;
            binWidth = 0.3;
          }
        }
        else if ( pTy1py2/myy >= 0.8 && pTy1py2/myy < 1.1 ) {
          if ( pTy1my2/myy < 0.05 ) { bin = 1; binWidth = 0.05;}
          else if ( pTy1my2/myy < 0.10) { bin = 2; binWidth = 0.05; }
          else if ( pTy1my2/myy < 0.20) { bin = 3; binWidth = 0.10; }
          else if ( pTy1my2/myy < 0.80) { bin = 4; binWidth = 0.60; }
        } else if ( pTy1py2/myy >= 0.8 && pTy1py2/myy < 4.0) {
          if ( pTy1my2/myy < 0.3 ) { bin = 5; binWidth = 0.3; }
          else if ( pTy1my2/myy < 0.6 ) { bin = 6; binWidth = 0.3; }
          else if ( pTy1my2/myy < 4.0 ) { bin = 7; binWidth = 3.4; }
        }
        if (bin>-1) {
          fillHist("rel_DpT_y_y_vs_rel_sumpT_y_y", bin+1, binWidth);
        }
      }

      // Jet multiplicity bin
      _h["N_j_30"]->fill( njets30>3 ? 4: njets30 + 1 );
      if (njets30central<1 || nlep10>0) 	_h["catXS_nbjet"]->fill(1.);
      else if (nbjets==0) 	_h["catXS_nbjet"]->fill(2.);
      else if (nbjets>=1) 	_h["catXS_nbjet"]->fill(3.);


      // Jet variables
      bool isVBF(false);
      FourMomentum j1(0.,0.,0.,0.);
      double pT_j1_30(0.);
      double maxtau (0.);
      double sumtau (0.);
      double pT_yyj_30 (0.);
      double m_yyj_30 (0.);
      double dphi_jj_signed(-99.);
      double ht = sum(jets30, Kin::pT, 0.0);
      _h["HT_30"]->fill(ht/GeV);
      if (njets30 > 0) {
        j1 = jets30[0].mom();
        pT_j1_30 = j1.pt();
        maxtau = max_tau_jet(yy, jets30);
        sumtau = sum_tau_jet(yy, jets30);
        pT_yyj_30 = (y1+y2+j1).pt();
        m_yyj_30 = (y1+y2+j1).mass();

        _h["pT_j1_30"]->fill(pT_j1_30 / GeV);
        _h["maxTau_yyj_30"]->fill(maxtau / GeV);
        _h["sumTau_yyj_30"]->fill(sumtau / GeV);
        _h["pT_yyj_30"]->fill(pT_yyj_30 / GeV);
        _h["m_yyj_30"]->fill(m_yyj_30 / GeV);

        if (pT_j1_30 < 40*GeV) _h["pT_yy_JV_40"]->fill(pT_yy / GeV);
        if (pT_j1_30 < 50*GeV) _h["pT_yy_JV_50"]->fill(pT_yy / GeV);
        if (pT_j1_30 < 60*GeV) _h["pT_yy_JV_60"]->fill(pT_yy / GeV);
      }
      else {
        _h["pT_j1_30"]->fill(15.) ;
        _h["pT_yy_JV_30"]->fill(pT_yy / GeV);
        _h["pT_yy_JV_40"]->fill(pT_yy / GeV);
        _h["pT_yy_JV_50"]->fill(pT_yy / GeV);
        _h["pT_yy_JV_60"]->fill(pT_yy / GeV);
      }

      // Dijet variables
      if ( njets30 > 1 ) {
        FourMomentum j2 = jets30[1].momentum();
        FourMomentum dijet = j1 + j2;
        FourMomentum yyjj = yy + dijet;
        double mjj = dijet.mass();
        dphi_jj_signed = j1.rap() > j2.rap()? deltaPhi(j1, j2, true) : deltaPhi(j2, j1, true);
        double dyjj = deltaRap(j1, j2);
        double dphiyyjj = deltaPhi(yy, dijet, false);
        double abs_Zepp = fabs( (y1+y2).eta() - 0.5*(j1.eta() + j2.eta()));

        _h["m_jj_30"]->fill(mjj / GeV);
        _h["Dphi_j_j_30_signed"]->fill(dphi_jj_signed);
        _h["pT_yyjj_30"]->fill(yyjj.pt() / GeV);
        _h["Dphi_yy_jj_30"]->fill(M_PI - dphiyyjj); // pi - dphi

        isVBF = (mjj/GeV>600. && dyjj>3.5);
        if (isVBF) {
          _h["VBF_pT_j1_30"]->fill(pT_j1_30 / GeV);
          _h["VBF_Dphi_j_j_30_signed"]->fill(dphi_jj_signed);
          _h["VBF_pT_yyjj_30"]->fill(yyjj.pt() / GeV);
          _h["VBF_abs_Zepp"]->fill(abs_Zepp);
        }
      }

      // 2d distributions

      // pT(yy) in bins pT(yyj)
      if (njets30==0) {
        if (pT_yy/GeV < 350) fillHist("pT_yy_vs_pT_yyj", 1, 350.);
      }
      else {
        if (pT_yyj_30/GeV < 30) {
          if (pT_yy/GeV < 100) fillHist("pT_yy_vs_pT_yyj", 2, 100.);
          else if (pT_yy/GeV < 350) fillHist("pT_yy_vs_pT_yyj", 3, 250.);
        }
        else if (pT_yyj_30/GeV < 60) {
          if (pT_yy/GeV < 45) fillHist("pT_yy_vs_pT_yyj", 4, 45.);
          else if (pT_yy/GeV < 120) fillHist("pT_yy_vs_pT_yyj", 5, 120.-45.);
          else if (pT_yy/GeV < 350) fillHist("pT_yy_vs_pT_yyj", 6, 350.-120.);
        }
        else if (pT_yyj_30/GeV < 350) {
          if (pT_yy/GeV < 80) fillHist("pT_yy_vs_pT_yyj", 7, 80.);
          else if (pT_yy/GeV < 250) fillHist("pT_yy_vs_pT_yyj", 8, 250.-80.);
          else if (pT_yy/GeV < 450) fillHist("pT_yy_vs_pT_yyj", 9, 450.-250.);
        }
      }

      // pT(yy) in bins of max(tau_C^j)
      if ( njets30 == 0 ) {
        if (pT_yy/GeV < 350.)
          fillHist("pT_yy_vs_maxTau_yyj", 1., 350.);
      }
      else {
        if ( maxtau/GeV < 15. ) {

          if (pT_yy/GeV < 100.)
            fillHist("pT_yy_vs_maxTau_yyj", 2., 100.);

          else if (pT_yy/GeV < 350.)
            fillHist("pT_yy_vs_maxTau_yyj", 3., 350.-100.);

        }
      	else if (maxtau/GeV < 25.) {

        if (pT_yy/GeV < 120.)
          fillHist("pT_yy_vs_maxTau_yyj", 4., 120.);

        else if (pT_yy/GeV < 350.)
          fillHist("pT_yy_vs_maxTau_yyj", 5., 350.-120.);

        }
      	else if (maxtau/GeV < 40.) {

        if (pT_yy/GeV < 200.)
          fillHist("pT_yy_vs_maxTau_yyj", 6., 200.);

        else if (pT_yy/GeV < 350.)
          fillHist("pT_yy_vs_maxTau_yyj", 7., 350.-200.);

        }
      	else if (maxtau/GeV < 400.) {

        if (pT_yy/GeV < 250.)
          fillHist("pT_yy_vs_maxTau_yyj", 8., 250. );

        else if (pT_yy/GeV < 650.)
          fillHist("pT_yy_vs_maxTau_yyj", 9., 650.-250.);

        }
      }

      // VBF pT_j1_30 in bins of Dphi_j_j_30_signed
      if (isVBF) {
        if (dphi_jj_signed < 0.) {

          if (pT_j1_30/GeV < 120.)
            fillHist("VBF_pT_j1_30_vs_Dphi_j_j_30_signed", 1., 120.-30.);

          else if (pT_j1_30/GeV < 500.)
            fillHist("VBF_pT_j1_30_vs_Dphi_j_j_30_signed", 2., 500.-120.);

        }
        else {

          if (pT_j1_30/GeV < 120.)
            fillHist("VBF_pT_j1_30_vs_Dphi_j_j_30_signed", 3., 120.-30.);

          else if (pT_j1_30/GeV < 500.)
            fillHist("VBF_pT_j1_30_vs_Dphi_j_j_30_signed", 4., 500.-120.);

        }
      }


      // fiducial regions
      // inclusive
      _h["fid_regions"]->fill( 1.0 ) ;
      // VBF
      if ( isVBF ) _h["fid_regions"]->fill( 2.0 );
      // lep
      if ( nlep15>0 ) _h["fid_regions"]->fill( 3.0 ) ;
      // MET
      if ( met>=80*GeV && pT_yy>80*GeV ) _h["fid_regions"]->fill( 4.0 ) ;
      // ttH
      const bool ttH_lep = njets30 > 2 && nlep15 >= 1 && nbjets > 0;
      const bool ttH_had = njets30 > 3 && nlep15 == 0 && nbjets > 0;
      if ( ttH_lep || ttH_had ) _h["fid_regions"]->fill( 5.0 );
    }


    void finalize () {
      // Scale histograms from nEvents (sumW) to Xsection
      const double xs = crossSectionPerEvent() / femtobarn;
      scale(_h, xs);
    }

    double tau_jet(const FourMomentum& Higgs, const Jet& jet ) const {
      const double mTj = sqrt( sqr(jet.pT()) + sqr(jet.mass()) );
      return mTj/(2.0*cosh( jet.rap() - Higgs.rap() ) ) ;
    }

    double max_tau_jet ( const FourMomentum& Higgs , const Jets& jets ) const {
      double max_tj ( 0 ) ;
      for (const Jet& jet : jets) {
        double tau_j ( tau_jet (Higgs , jet) ) ;
        if ( tau_j < max_tj ) continue ;
        max_tj = tau_j ;
      }
      return max_tj ;
    }

    double sum_tau_jet ( const FourMomentum& Higgs , const Jets& jets) const {
      double sum_tj ( 0 ) ;
      double temp ( -99 ) ;
      for (const Jet& jet : jets){
        temp =  tau_jet(Higgs, jet);
        //sum_tj += temp;
        sum_tj += temp > 8*GeV ? temp : 0.0;
      }
      return sum_tj ;
    }

    void fillHist(const string& name, const double value, const double binWidth=-1.0) {
      // scale weight by 1/binWidth to make differential xsections distributions
      // hack needed for 2D hists which are shown as 1D histograms
      if (binWidth) _h[name]->fill(value , 1. / binWidth);
      else {
        // handle overflow/underflow
        if (value > _h[name]->xMin() && value < _h[name]->xMax()) {
          _h[name]->fill(value , 1. / _h[name]->binAt(value).xWidth());
        }
        else  _h[name]->fill(value);
      }
    }

  private:

    map<string,Histo1DPtr> _h;

    //@}
  } ;

  // This is a required hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2022_I2023464);
}
