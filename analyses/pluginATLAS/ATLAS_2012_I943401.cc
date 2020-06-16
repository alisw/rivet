// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  class ATLAS_2012_I943401 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor

    ATLAS_2012_I943401()
      : Analysis("ATLAS_2012_I943401")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 20*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // jet finder
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3 && Cuts::pT > 0.5*GeV), "cfs");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.5), "vfs");

      // book histograms

      // counts in signal regions
      book(_count_OS_SR1 ,"count_OS_SR1", 1, 0., 1.);
      book(_count_OS_SR2 ,"count_OS_SR2", 1, 0., 1.);
      book(_count_OS_SR3 ,"count_OS_SR3", 1, 0., 1.);
      book(_count_SS_SR1 ,"count_SS_SR1", 1, 0., 1.);
      book(_count_SS_SR2 ,"count_SS_SR2", 1, 0., 1.);
      book(_count_FS_SR1 ,"count_FS_SR1", 1, 0., 1.);
      book(_count_FS_SR2 ,"count_FS_SR2", 1, 0., 1.);
      book(_count_FS_SR3 ,"count_FS_SR3", 1, 0., 1.);

      // histograms from paper

      book(_hist_mll_SS_D         , 1,1,1);
      book(_hist_mll_SS_B         , 1,1,2);
      book(_hist_eTmiss_SS_D      , 2,1,1);
      book(_hist_eTmiss_SS_B      , 2,1,2);
      book(_hist_mll_SS_2Jet_D    , 3,1,1);
      book(_hist_mll_SS_2Jet_B    , 3,1,2);
      book(_hist_njet_SS_D        , 5,1,1);
      book(_hist_njet_SS_B        , 5,1,2);
      book(_hist_pT_j1_SS_D       , 6,1,1);
      book(_hist_pT_j1_SS_B       , 6,1,2);
      book(_hist_pT_j2_SS_D       , 7,1,1);
      book(_hist_pT_j2_SS_B       , 7,1,2);
      book(_hist_pT_l1_SS_D       , 8,1,1);
      book(_hist_pT_l1_SS_B       , 8,1,2);
      book(_hist_pT_l2_SS_D       , 9,1,1);
      book(_hist_pT_l2_SS_B       , 9,1,2);
      book(_hist_mll_OS_D         ,10,1,1);
      book(_hist_mll_OS_B         ,10,1,2);
      book(_hist_eTmiss_OS_D      ,11,1,1);
      book(_hist_eTmiss_OS_B      ,11,1,2);
      book(_hist_eTmiss_3Jet_OS_D ,12,1,1);
      book(_hist_eTmiss_3Jet_OS_B ,12,1,2);
      book(_hist_eTmiss_4Jet_OS_D ,13,1,1);
      book(_hist_eTmiss_4Jet_OS_B ,13,1,2);
      book(_hist_njet_OS_D        ,14,1,1);
      book(_hist_njet_OS_B        ,14,1,2);
      book(_hist_pT_j1_OS_D       ,15,1,1);
      book(_hist_pT_j1_OS_B       ,15,1,2);
      book(_hist_pT_j2_OS_D       ,16,1,1);
      book(_hist_pT_j2_OS_B       ,16,1,2);
      book(_hist_pT_l1_OS_D       ,17,1,1);
      book(_hist_pT_l1_OS_B       ,17,1,2);
      book(_hist_pT_l2_OS_D       ,18,1,1);
      book(_hist_pT_l2_OS_B       ,18,1,2);
      //????
      //   <dataPointSet name="d04-x01-y01" dimension="2" path="/REF/ATLAS_2011_I943401" title="EVENTS/10 GEV" >
      //   <dataPointSet name="d04-x01-y02" dimension="2" path="/REF/ATLAS_2011_I943401" title="EVENTS/10 GEV" >
    }

    /// Perform the event analysis
    void analyze(const Event& event) {
      // get the jet candidates
      Jets cand_jets;
      for (const Jet& jet :
        apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.eta() ) < 2.8 ) {
          cand_jets.push_back(jet);
        }
      }

      // electron candidates
      Particles cand_e =
        apply<IdentifiedFinalState>(event, "elecs").particlesByPt();

      // Discard jets that overlap with electrons
      Jets recon_jets;
      for ( const Jet& jet : cand_jets ) {
        bool away_from_e = true;
          for ( const Particle& e : cand_e ) {
            if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
              away_from_e = false;
              break;
            }
          }
          if ( away_from_e ) recon_jets.push_back( jet );
      }
      // get the charged tracks for isolation
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();

      // Reconstructed electrons
      Particles recon_e;
      for ( const Particle& e : cand_e ) {
        // check not near a jet
        bool e_near_jet = false;
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
            e_near_jet = true;
            break;
          }
        }
        if ( e_near_jet ) continue;
        // check the isolation
        double pTinCone = -e.pT();
        for ( const Particle& track : chg_tracks ) {
          if ( deltaR(e.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 0.1*e.perp() )
          recon_e.push_back(e);
      }

      // Reconstructed Muons
      Particles recon_mu;
      Particles cand_mu =
        apply<IdentifiedFinalState>(event,"muons").particlesByPt();
      for ( const Particle& mu : cand_mu ) {
        // check not near a jet
        bool mu_near_jet = false;
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            mu_near_jet = true;
            break;
          }
        }
        if ( mu_near_jet ) continue;
        // isolation
        double pTinCone = -mu.pT();
        for ( const Particle& track : chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 1.8*GeV )
          recon_mu.push_back(mu);
      }

      // pTmiss
      Particles vfs_particles
        = apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle& p : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // ATLAS calo problem
      if(rand()/static_cast<double>(RAND_MAX)<=0.42) {
        for ( const Particle& e : recon_e ) {
          double eta = e.eta();
          double phi = e.azimuthalAngle(MINUSPI_PLUSPI);
          if (inRange(eta, -0.1, 1.5) && inRange(phi, -0.9, -0.5)) vetoEvent;
        }
        for ( const Jet& jet : recon_jets ) {
          double eta = jet.rapidity();
          double phi = jet.azimuthalAngle(MINUSPI_PLUSPI);
          if (jet.pT() > 40*GeV && inRange(eta, -0.1, 1.5) && inRange(phi, -0.9, -0.5)) vetoEvent;
        }
      }

      // Exactly two leptons for each event
      if ( recon_mu.size() + recon_e.size() != 2) vetoEvent;
      // two electrons highest pT > 25
      Particles recon_leptons;
      if (recon_e.size()==2 && recon_e[0].pT()>25*GeV) {
        recon_leptons = recon_e;
      }
      // two muons highest pT > 20
      else if (recon_mu.size()==2 && recon_mu[0].pT() > 20*GeV) {
        recon_leptons = recon_mu;
      } else if (recon_e.size()==1 && recon_mu.size()==1 &&
                 (recon_e[0].pT() > 25*GeV || recon_mu[0].pT() > 20*GeV )) {
        if (recon_mu[0].pT() < recon_e[0].pT()) {
          recon_leptons.push_back(recon_e [0]);
          recon_leptons.push_back(recon_mu[0]);
        } else {
          recon_leptons.push_back(recon_mu[0]);
          recon_leptons.push_back(recon_e [0]);
        }
      }
      // fails trigger
      else vetoEvent;

      double mll = (recon_leptons[0].momentum()+recon_leptons[1].momentum()).mass();
      // lepton pair mass > 12.
      if (mll < 12*GeV) vetoEvent;

      // same sign or opposite sign event
      int sign = recon_leptons[0].pid()*recon_leptons[1].pid();

      // same sign leptons
      if(sign>0) {
        _hist_mll_SS_D   ->fill(mll   );
        _hist_mll_SS_B   ->fill(mll   );
        _hist_eTmiss_SS_D->fill(eTmiss);
        _hist_eTmiss_SS_B->fill(eTmiss);
        if(recon_jets.size()>=2) {
          _hist_mll_SS_2Jet_D   ->fill(mll   );
          _hist_mll_SS_2Jet_B   ->fill(mll   );
        }
        _hist_njet_SS_D ->fill(recon_jets.size());
        _hist_njet_SS_B ->fill(recon_jets.size());
        if(!recon_jets.empty()) {
          _hist_pT_j1_SS_D->fill(recon_jets[0].perp());
          _hist_pT_j1_SS_B->fill(recon_jets[0].perp());
        }
        if(recon_jets.size()>2) {
          _hist_pT_j2_SS_D->fill(recon_jets[1].perp());
          _hist_pT_j2_SS_B->fill(recon_jets[1].perp());
        }
        _hist_pT_l1_SS_D->fill(recon_leptons[0].perp());
        _hist_pT_l1_SS_B->fill(recon_leptons[0].perp());
        _hist_pT_l2_SS_D->fill(recon_leptons[1].perp());
        _hist_pT_l2_SS_B->fill(recon_leptons[1].perp());
        // SS-SR1
        if(eTmiss>100.) {
          _count_SS_SR1->fill(0.5);
        }
        // SS-SR2
        if(eTmiss>80. && recon_jets.size()>=2 &&
           recon_jets[1].perp()>50.) {
          _count_SS_SR2->fill(0.5);
        }
      }
      // opposite sign
      else {
        _hist_mll_OS_D->fill(mll   );
        _hist_mll_OS_B->fill(mll   );
        _hist_eTmiss_OS_D->fill(eTmiss);
        _hist_eTmiss_OS_B->fill(eTmiss);
        if(recon_jets.size()>=3){
          _hist_eTmiss_3Jet_OS_D->fill(eTmiss);
          _hist_eTmiss_3Jet_OS_B->fill(eTmiss);
        }
        if(recon_jets.size()>=4){
          _hist_eTmiss_4Jet_OS_D->fill(eTmiss);
          _hist_eTmiss_4Jet_OS_B->fill(eTmiss);
        }
        _hist_njet_OS_D->fill(recon_jets.size());
        _hist_njet_OS_B->fill(recon_jets.size());
        if(!recon_jets.empty()) {
          _hist_pT_j1_OS_D->fill(recon_jets[0].perp());
          _hist_pT_j1_OS_B->fill(recon_jets[0].perp());
        }
        if(recon_jets.size()>2) {
          _hist_pT_j2_OS_D->fill(recon_jets[1].perp());
          _hist_pT_j2_OS_B->fill(recon_jets[1].perp());
        }
        _hist_pT_l1_OS_D->fill(recon_leptons[0].perp());
        _hist_pT_l1_OS_B->fill(recon_leptons[0].perp());
        _hist_pT_l2_OS_D->fill(recon_leptons[1].perp());
        _hist_pT_l2_OS_B->fill(recon_leptons[1].perp());
        // different signal regions
        // OS-SR1
        if(eTmiss>250.) {
          _count_OS_SR1->fill(0.5);
        }
        // OS-SR2
        if(eTmiss>220. && recon_jets.size()>=3 &&
           recon_jets[0].perp()>80. &&
           recon_jets[2].perp()>40.) {
          _count_OS_SR2->fill(0.5);
        }
        // OS-SR3
        if(eTmiss>100. && recon_jets.size()>=4 &&
           recon_jets[0].perp()>100. &&
           recon_jets[3].perp()>70.) {
          _count_OS_SR3->fill(0.5);
        }
        // same flavour analysis
        static const double beta   = 0.75;
        static const double tau_e  = 0.96;
        static const double tau_mu = 0.816;
        double fs_weight = 1.0;
        if (recon_leptons[0].abspid() == PID::ELECTRON && recon_leptons[1].abspid() == PID::ELECTRON) {
          fs_weight /= beta*(1.-sqr(1.-tau_e));
        } else if (recon_leptons[0].abspid() == PID::MUON && recon_leptons[1].abspid()==PID::MUON) {
          fs_weight *= beta/(1.-sqr(1.-tau_mu));
        } else {
          fs_weight /= -(1.-(1.-tau_e)*(1.-tau_mu));
        }
        // FS-SR1
        if(eTmiss>80.&& (mll<80.||mll>100.)) {
          _count_FS_SR1->fill(0.5,fs_weight);
        }
        // FS-SR2
        if(eTmiss>80.&&recon_jets.size()>=2) {
          _count_FS_SR2->fill(0.5,fs_weight);
        }
        // FS-SR3
        if(eTmiss>250.) {
          _count_FS_SR3->fill(0.5,fs_weight);
        }
      }
    }

    //@}


    void finalize() {

      double norm = crossSection()/femtobarn*1.04/sumOfWeights();
      // event counts
      scale(_count_OS_SR1,norm);
      scale(_count_OS_SR2,norm);
      scale(_count_OS_SR3,norm);
      scale(_count_SS_SR1,norm);
      scale(_count_SS_SR2,norm);
      scale(_count_FS_SR1,norm);
      scale(_count_FS_SR2,norm);
      scale(_count_FS_SR3,norm);
      // histograms
      scale(_hist_mll_SS_D     ,norm*20.);
      scale(_hist_mll_SS_B     ,norm*20.);
      scale(_hist_eTmiss_SS_D  ,norm*20.);
      scale(_hist_eTmiss_SS_B  ,norm*20.);
      scale(_hist_mll_SS_2Jet_D,norm*50.);
      scale(_hist_mll_SS_2Jet_B,norm*50.);
      scale(_hist_njet_SS_D    ,norm    );
      scale(_hist_njet_SS_B    ,norm    );
      scale(_hist_pT_j1_SS_D   ,norm*20.);
      scale(_hist_pT_j1_SS_B   ,norm*20.);
      scale(_hist_pT_j2_SS_D   ,norm*20.);
      scale(_hist_pT_j2_SS_B   ,norm*20.);
      scale(_hist_pT_l1_SS_D   ,norm*5. );
      scale(_hist_pT_l1_SS_B   ,norm*5. );
      scale(_hist_pT_l2_SS_D   ,norm*5. );
      scale(_hist_pT_l2_SS_B   ,norm*5. );

      scale(_hist_mll_OS_D        ,norm*10.);
      scale(_hist_mll_OS_B        ,norm*10.);
      scale(_hist_eTmiss_OS_D     ,norm*10.);
      scale(_hist_eTmiss_OS_B     ,norm*10.);
      scale(_hist_eTmiss_3Jet_OS_D,norm*10.);
      scale(_hist_eTmiss_3Jet_OS_B,norm*10.);
      scale(_hist_eTmiss_4Jet_OS_D,norm*10.);
      scale(_hist_eTmiss_4Jet_OS_B,norm*10.);
      scale(_hist_njet_OS_D       ,norm    );
      scale(_hist_njet_OS_B       ,norm    );
      scale(_hist_pT_j1_OS_D      ,norm*20.);
      scale(_hist_pT_j1_OS_B      ,norm*20.);
      scale(_hist_pT_j2_OS_D      ,norm*20.);
      scale(_hist_pT_j2_OS_B      ,norm*20.);
      scale(_hist_pT_l1_OS_D      ,norm*20.);
      scale(_hist_pT_l1_OS_B      ,norm*20.);
      scale(_hist_pT_l2_OS_D      ,norm*20.);
      scale(_hist_pT_l2_OS_B      ,norm*20.);
    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _count_OS_SR1;
    Histo1DPtr _count_OS_SR2;
    Histo1DPtr _count_OS_SR3;
    Histo1DPtr _count_SS_SR1;
    Histo1DPtr _count_SS_SR2;
    Histo1DPtr _count_FS_SR1;
    Histo1DPtr _count_FS_SR2;
    Histo1DPtr _count_FS_SR3;

    Histo1DPtr _hist_mll_SS_D;
    Histo1DPtr _hist_mll_SS_B;
    Histo1DPtr _hist_eTmiss_SS_D;
    Histo1DPtr _hist_eTmiss_SS_B;
    Histo1DPtr _hist_mll_SS_2Jet_D;
    Histo1DPtr _hist_mll_SS_2Jet_B;
    Histo1DPtr _hist_njet_SS_D;
    Histo1DPtr _hist_njet_SS_B;
    Histo1DPtr _hist_pT_j1_SS_D;
    Histo1DPtr _hist_pT_j1_SS_B;
    Histo1DPtr _hist_pT_j2_SS_D;
    Histo1DPtr _hist_pT_j2_SS_B;
    Histo1DPtr _hist_pT_l1_SS_D;
    Histo1DPtr _hist_pT_l1_SS_B;
    Histo1DPtr _hist_pT_l2_SS_D;
    Histo1DPtr _hist_pT_l2_SS_B;

    Histo1DPtr _hist_mll_OS_D;
    Histo1DPtr _hist_mll_OS_B;
    Histo1DPtr _hist_eTmiss_OS_D;
    Histo1DPtr _hist_eTmiss_OS_B;
    Histo1DPtr _hist_eTmiss_3Jet_OS_D;
    Histo1DPtr _hist_eTmiss_3Jet_OS_B;
    Histo1DPtr _hist_eTmiss_4Jet_OS_D;
    Histo1DPtr _hist_eTmiss_4Jet_OS_B;
    Histo1DPtr _hist_njet_OS_D ;
    Histo1DPtr _hist_njet_OS_B ;
    Histo1DPtr _hist_pT_j1_OS_D;
    Histo1DPtr _hist_pT_j1_OS_B;
    Histo1DPtr _hist_pT_j2_OS_D;
    Histo1DPtr _hist_pT_j2_OS_B;
    Histo1DPtr _hist_pT_l1_OS_D;
    Histo1DPtr _hist_pT_l1_OS_B;
    Histo1DPtr _hist_pT_l2_OS_D;
    Histo1DPtr _hist_pT_l2_OS_B;
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I943401);

}
