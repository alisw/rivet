// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2011_CONF_2011_090 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor

    ATLAS_2011_CONF_2011_090()
      : Analysis("ATLAS_2011_CONF_2011_090")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialize projections before the run
    void init() {

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT >= 20*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // veto region electrons (from 2010 arXiv:1102.2357v2)
      Cut vetocut = Cuts::absetaIn(1.37, 1.52);
      IdentifiedFinalState veto_elecs(vetocut && Cuts::pT > 10*GeV);
      veto_elecs.acceptIdPair(PID::ELECTRON);
      declare(veto_elecs, "veto_elecs");

      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // Jet finder
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3.0 && Cuts::pT > 0.5*GeV), "cfs");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9), "vfs");

      /// Book histograms
      book(_count_mu_channel ,"count_muon_channel", 1, 0., 1.);
      book(_count_e_channel ,"count_electron_channel", 1, 0., 1.);
      book(_hist_eTmiss_e ,"Et_miss_e", 50, 0., 500.);
      book(_hist_eTmiss_mu ,"Et_miss_mu", 50, 0., 500.);
      book(_hist_m_eff_e ,"m_eff_e", 60, 0., 1500.);
      book(_hist_m_eff_mu ,"m_eff_mu", 60, 0., 1500.);
      book(_hist_m_eff_e_final ,"m_eff_e_final", 15, 0., 1500.);
      book(_hist_m_eff_mu_final ,"m_eff_mu_final", 15, 0., 1500.);
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = 1.0;

      Particles veto_e = apply<IdentifiedFinalState>(event, "veto_elecs").particles();
      if ( ! veto_e.empty() ) {
        MSG_DEBUG("electrons in veto region");
        vetoEvent;
      }

      Jets cand_jets = apply<FastJets>(event, "AntiKtJets04").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.8);

      Particles candtemp_e = apply<IdentifiedFinalState>(event, "elecs").particlesByPt();
      Particles candtemp_mu = apply<IdentifiedFinalState>(event,"muons").particlesByPt();
      Particles chg_tracks = apply<ChargedFinalState>(event, "cfs").particles();
      Particles cand_mu;
      Particles cand_e;

      // pTcone around muon track
      for ( const Particle& mu : candtemp_mu ) {
        double pTinCone = -mu.pT();
        for ( const Particle& track : chg_tracks ) {
          if ( deltaR(mu, track) < 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 1.8*GeV )
          cand_mu.push_back(mu);
      }


      // pTcone around electron
      for ( const Particle& e : candtemp_e ) {
        double pTinCone = -e.pT();
        for ( const Particle& track : chg_tracks ) {
          if ( deltaR(e, track) < 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 0.10 * e.pT() )
          cand_e.push_back(e);
      }


      // discard jets that overlap with electrons
      Jets cand_jets_2;
      for ( const Jet& jet : cand_jets ) {
        bool away_from_e = true;
        for ( const Particle & e : cand_e ) {
          if ( deltaR(e, jet) <= 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e )
          cand_jets_2.push_back( jet );
      }

      // only consider leptons far from jet
      Particles recon_e, recon_mu;
      for ( const Particle & e : cand_e ) {
        bool e_near_jet = false;
        for ( const Jet& jet : cand_jets_2 ) {
          if (inRange(deltaR(e, jet), 0.2, 0.4)) e_near_jet = true;
        }
        if ( ! e_near_jet ) recon_e.push_back( e );
      }

      for ( const Particle & mu : cand_mu ) {
        bool mu_near_jet = false;
        for ( const Jet& jet : cand_jets_2 ) {
          if ( deltaR(mu, jet) < 0.4 ) mu_near_jet = true;
        }
        if ( ! mu_near_jet ) recon_mu.push_back( mu );
      }

      // pTmiss
      Particles vfs_particles
        = apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle & p : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();


      // final jet filter
      Jets recon_jets;
      for ( const Jet& jet : cand_jets_2 ) {
        recon_jets.push_back( jet );
      }




      // ==================== observables ====================


      // Njets

      int Njets = 0;
      double pTmiss_phi = pTmiss.phi();
      for ( const Jet& jet : recon_jets ) {
        if ( jet.abseta() < 2.8 )
          Njets+=1;
      }
      if ( Njets < 3 ) {
        MSG_DEBUG("Only " << Njets << " jets w/ eta<2.8 left");
        vetoEvent;
      }

      if ( recon_jets[0].pT() <= 60.0 * GeV ) {
        MSG_DEBUG("No hard leading jet in " << recon_jets.size() << " jets");
        vetoEvent;
      }
      for ( int i = 1; i < 3; ++i ) {
        if ( recon_jets[i].pT() <= 25*GeV ) {
          vetoEvent;
        }
      }

      for ( int i = 0; i < 3; ++i ) {
        double dPhi = deltaPhi( pTmiss_phi, recon_jets[i].phi() );
        if ( dPhi <= 0.2 ) {
          MSG_DEBUG("dPhi too small");
          vetoEvent;
          break;
        }
      }


      Particles lepton;
      if ( recon_mu.empty() && recon_e.empty() ) {
        MSG_DEBUG("No leptons");
        vetoEvent;
      }
      else {
        for ( const Particle & mu : recon_mu )
          lepton.push_back(mu);
        for ( const Particle & e : recon_e )
          lepton.push_back(e);
      }

      std::sort(lepton.begin(), lepton.end(), cmpMomByPt);

      double e_id = 11;
      double mu_id = 13;

      // one hard leading lepton cut
      if ( lepton[0].abspid() == e_id &&
           lepton[0].pT() <= 25*GeV ) {
        vetoEvent;
      }
      else if ( lepton[0].abspid() == mu_id &&
                lepton[0].pT() <= 20*GeV ) {
        vetoEvent;
      }

      // exactly one hard leading lepton cut
      if(lepton.size()>1) {
        if ( lepton[1].abspid() == e_id &&
             lepton[1].pT() > 20*GeV ) {
          vetoEvent;
        }
        else if ( lepton[1].abspid() == mu_id &&
                  lepton[1].pT() > 10*GeV ) {
          vetoEvent;
        }
      }


      // ==================== FILL ====================


      FourMomentum pT_l = lepton[0].momentum();


      double dPhi = deltaPhi( pT_l.phi(), pTmiss_phi);
      double mT = sqrt( 2 * pT_l.pT() * eTmiss * (1 - cos(dPhi)) );


      // effective mass
      double m_eff = eTmiss + pT_l.pT()
        + recon_jets[0].pT()
        + recon_jets[1].pT()
        + recon_jets[2].pT();


      // Electron channel signal region

      if (  lepton[0].abspid() == e_id ) {

        _hist_eTmiss_e->fill(eTmiss, weight);
        _hist_m_eff_e->fill(m_eff, weight);

        if ( mT > 100*GeV && eTmiss > 125*GeV ) {
          _hist_m_eff_e_final->fill(m_eff, weight);
          if ( m_eff > 500*GeV && eTmiss > 0.25*m_eff ) {
            _count_e_channel->fill(0.5,weight);
          }
        }
      }

      // Muon channel signal region

      else if (  lepton[0].abspid() == mu_id ) {

        _hist_eTmiss_mu->fill(eTmiss, weight);
        _hist_m_eff_mu->fill(m_eff, weight);

        if ( mT > 100*GeV && eTmiss > 125*GeV ) {
          _hist_m_eff_mu_final->fill(m_eff, weight);
          if ( m_eff > 500*GeV && eTmiss > 0.25*m_eff ) {
            _count_mu_channel->fill(0.5,weight);
          }
        }

      }


    }

    //@}


    void finalize() {
      scale( _hist_eTmiss_e      , 10.  * 165. * crossSection()/picobarn/sumOfWeights() );
      scale( _hist_eTmiss_mu     , 10.  * 165. * crossSection()/picobarn/sumOfWeights() );
      scale( _hist_m_eff_e       , 25.  * 165. * crossSection()/picobarn/sumOfWeights() );
      scale( _hist_m_eff_mu      , 25.  * 165. * crossSection()/picobarn/sumOfWeights() );
      scale( _hist_m_eff_e_final , 100. * 165. * crossSection()/picobarn/sumOfWeights() );
      scale( _hist_m_eff_mu_final, 100. * 165. * crossSection()/picobarn/sumOfWeights() );
    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _count_e_channel;
    Histo1DPtr _count_mu_channel;

    Histo1DPtr _hist_eTmiss_e;
    Histo1DPtr _hist_eTmiss_mu;

    Histo1DPtr _hist_m_eff_e;
    Histo1DPtr _hist_m_eff_mu;
    Histo1DPtr _hist_m_eff_e_final;
    Histo1DPtr _hist_m_eff_mu_final;


    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_CONF_2011_090);

}
