// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {



  /// @author Peter Richardson
  class ATLAS_2012_I1125961 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1125961()
      : Analysis("ATLAS_2012_I1125961")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 20*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // Projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // Jet finder
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // All tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3.0), "cfs");

      // Used for pTmiss (N.B. the real 'vfs' extends beyond 4.5 to |eta| = 4.9)
      declare(VisibleFinalState(Cuts::abseta < 4.5), "vfs");

      // Book histograms
      book(_count_A_tight   ,"count_A_tight"   , 1, 0., 1.);
      book(_count_A_medium  ,"count_A_medium"  , 1, 0., 1.);
      book(_count_Ap_medium ,"count_Ap_medium" , 1, 0., 1.);
      book(_count_B_tight   ,"count_B_tight"   , 1, 0., 1.);
      book(_count_C_tight   ,"count_C_tight"   , 1, 0., 1.);
      book(_count_C_medium  ,"count_C_medium"  , 1, 0., 1.);
      book(_count_C_loose   ,"count_C_loose"   , 1, 0., 1.);
      book(_count_D_tight   ,"count_D_tight"   , 1, 0., 1.);
      book(_count_E_tight   ,"count_E_tight"   , 1, 0., 1.);
      book(_count_E_medium  ,"count_E_medium"  , 1, 0., 1.);
      book(_count_E_loose   ,"count_E_loose"   , 1, 0., 1.);

      book(_hist_meff_A  ,"hist_m_eff_A" , 30, 0., 3000.);
      book(_hist_meff_Ap ,"hist_m_eff_Ap", 30, 0., 3000.);
      book(_hist_meff_B  ,"hist_m_eff_B" , 30, 0., 3000.);
      book(_hist_meff_C  ,"hist_m_eff_C" , 30, 0., 3000.);
      book(_hist_meff_D  ,"hist_m_eff_D" , 30, 0., 3000.);
      book(_hist_meff_E  ,"hist_m_eff_E" , 30, 0., 3000.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      Jets cand_jets;
      const Jets jets = apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV);
      for (const Jet& jet : jets) {
        if ( fabs( jet.eta() ) < 4.9 ) {
          cand_jets.push_back(jet);
        }
      }

      const Particles cand_e  = apply<IdentifiedFinalState>(event, "elecs").particlesByPt();

      // Muon isolation not mentioned in hep-exp 1109.6572 but assumed to still be applicable
      Particles cand_mu;
      const Particles chg_tracks = apply<ChargedFinalState>(event, "cfs").particles();
      const Particles muons = apply<IdentifiedFinalState>(event, "muons").particlesByPt();
      for (const Particle& mu : muons) {
        double pTinCone = -mu.pT();
        for (const Particle& track : chg_tracks) {
          if ( deltaR(mu.momentum(),track.momentum()) <= 0.2 ) {
            pTinCone += track.pT();
          }
        }
        if ( pTinCone < 1.8*GeV ) cand_mu.push_back(mu);
      }

      // Resolve jet-lepton overlap for jets with |eta| < 2.8
      Jets recon_jets;
      for ( const Jet& jet : cand_jets ) {
        if ( fabs( jet.eta() ) >= 2.8 ) continue;
        bool away_from_e = true;
        for ( const Particle & e : cand_e ) {
          if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e ) recon_jets.push_back( jet );
      }

      Particles recon_e, recon_mu;

      for ( const Particle & e : cand_e ) {
        bool away = true;
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away ) recon_e.push_back( e );
      }

      for ( const Particle & mu : cand_mu ) {
        bool away = true;
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away ) recon_mu.push_back( mu );
      }

      // pTmiss
      // Based on all candidate electrons, muons and jets, plus everything else with |eta| < 4.5
      // i.e. everything in our projection "vfs" plus the jets with |eta| > 4.5
      Particles vfs_particles = apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle & p : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      for ( const Jet& jet : cand_jets ) {
        if ( fabs( jet.eta() ) > 4.5 ) pTmiss -= jet.momentum();
      }
      double eTmiss = pTmiss.pT();

      // no electron pT> 20 or muons pT>10
      if ( !recon_mu.empty() || !recon_e.empty() ) {
        MSG_DEBUG("Charged leptons left after selection");
        vetoEvent;
      }

      if ( eTmiss <= 160 * GeV ) {
        MSG_DEBUG("Not enough eTmiss: " << eTmiss << " < 130");
        vetoEvent;
      }

      if ( recon_jets.size()<2 ||
           recon_jets[0].pT() <= 130.0 * GeV ||
           recon_jets[0].pT() <=  60.0 * GeV ) {
        MSG_DEBUG("No hard leading jet in " << recon_jets.size() << " jets");
        vetoEvent;
      }

      // ==================== observables ====================

      int Njets = 0;
      double min_dPhi_All = 999.999;
      double min_dPhi_2   = 999.999;
      double min_dPhi_3   = 999.999;
      double pTmiss_phi = pTmiss.phi();
      for ( const Jet& jet : recon_jets ) {
        if ( jet.pT() < 40 * GeV ) continue;
        if ( Njets < 2 ) {
          min_dPhi_2 = min( min_dPhi_2, deltaPhi( pTmiss_phi, jet.phi() ) );
        }
        if( Njets < 3) {
          min_dPhi_3 = min( min_dPhi_3, deltaPhi( pTmiss_phi, jet.phi() ) );
        }
        min_dPhi_All = min( min_dPhi_All, deltaPhi( pTmiss_phi, jet.phi() ) );
        ++Njets;
      }

      // inclusive meff
      double m_eff_inc =  eTmiss;
      for ( const Jet& jet : recon_jets ) {
        double perp =  jet.pT();
        if(perp>40.) m_eff_inc += perp;
      }

      // region A
      double m_eff_Nj = eTmiss + recon_jets[0].pT() + recon_jets[1].pT();
      if( min_dPhi_2 > 0.4 && eTmiss/m_eff_Nj > 0.3 ) {
        _hist_meff_A->fill(m_eff_inc,weight);
        if(m_eff_inc>1900.) _count_A_tight ->fill(0.5,weight);
        if(m_eff_inc>1400.) _count_A_medium->fill(0.5,weight);
      }

      // region A'
      if( min_dPhi_2 > 0.4 && eTmiss/m_eff_Nj > 0.4 ) {
        _hist_meff_Ap->fill(m_eff_inc,weight);
        if(m_eff_inc>1200.) _count_Ap_medium->fill(0.5,weight);
      }

      // for rest of regions 3 jets pT> 60 needed
      if(recon_jets.size()<3 || recon_jets[2].perp()<60.)
        vetoEvent;

      // region B
      m_eff_Nj +=  recon_jets[2].perp();
      if( min_dPhi_3 > 0.4 && eTmiss/m_eff_Nj > 0.25 ) {
        _hist_meff_B->fill(m_eff_inc,weight);
        if(m_eff_inc>1900.) _count_B_tight ->fill(0.5,weight);
      }

      // for rest of regions 4 jets pT> 60 needed
      if(recon_jets.size()<4 || recon_jets[3].perp()<60.)
        vetoEvent;

      // region C
      m_eff_Nj +=  recon_jets[3].perp();
      if( min_dPhi_3 > 0.4 && min_dPhi_All > 0.2 && eTmiss/m_eff_Nj > 0.25 ) {
        _hist_meff_C->fill(m_eff_inc,weight);
        if(m_eff_inc>1500.) _count_C_tight ->fill(0.5,weight);
        if(m_eff_inc>1200.) _count_C_medium->fill(0.5,weight);
        if(m_eff_inc> 900.) _count_C_loose ->fill(0.5,weight);
      }

      // for rest of regions 5 jets pT> 40 needed
      if(recon_jets.size()<5 || recon_jets[4].perp()<40.)
        vetoEvent;

      // region D
      m_eff_Nj +=  recon_jets[4].perp();
      if( min_dPhi_3 > 0.4 && min_dPhi_All > 0.2 && eTmiss/m_eff_Nj > 0.2 ) {
        _hist_meff_D->fill(m_eff_inc,weight);
        if(m_eff_inc>1500.) _count_D_tight ->fill(0.5,weight);
      }

      // for rest of regions 6 jets pT> 40 needed
      if(recon_jets.size()<6 || recon_jets[5].perp()<40.)
        vetoEvent;

      // region E
      m_eff_Nj +=  recon_jets[5].perp();
      if( min_dPhi_3 > 0.4 && min_dPhi_All > 0.2 && eTmiss/m_eff_Nj > 0.15 ) {
        _hist_meff_E->fill(m_eff_inc,weight);
        if(m_eff_inc>1400.) _count_E_tight ->fill(0.5,weight);
        if(m_eff_inc>1200.) _count_E_medium->fill(0.5,weight);
        if(m_eff_inc> 900.) _count_E_loose ->fill(0.5,weight);
      }
    }


    void finalize() {

      double norm = crossSection()/femtobarn*4.7/sumOfWeights();
      // these are number of events at 4.7fb^-1 per 100 GeV
      scale( _hist_meff_A , 100. * norm );
      scale( _hist_meff_Ap, 100. * norm );
      scale( _hist_meff_B , 100. * norm );
      scale( _hist_meff_C , 100. * norm );
      scale( _hist_meff_D , 100. * norm );
      scale( _hist_meff_E , 100. * norm );
      // these are number of events at 4.7fb^-1
      scale(_count_A_tight  ,norm);
      scale(_count_A_medium ,norm);
      scale(_count_Ap_medium,norm);
      scale(_count_B_tight  ,norm);
      scale(_count_C_tight  ,norm);
      scale(_count_C_medium ,norm);
      scale(_count_C_loose  ,norm);
      scale(_count_D_tight  ,norm);
      scale(_count_E_tight  ,norm);
      scale(_count_E_medium ,norm);
      scale(_count_E_loose  ,norm);
    }

    //@}


  private:

    Histo1DPtr _count_A_tight;
    Histo1DPtr _count_A_medium;
    Histo1DPtr _count_Ap_medium;
    Histo1DPtr _count_B_tight;
    Histo1DPtr _count_C_tight;
    Histo1DPtr _count_C_medium;
    Histo1DPtr _count_C_loose;
    Histo1DPtr _count_D_tight;
    Histo1DPtr _count_E_tight;
    Histo1DPtr _count_E_medium;
    Histo1DPtr _count_E_loose;

    Histo1DPtr _hist_meff_A ;
    Histo1DPtr _hist_meff_Ap;
    Histo1DPtr _hist_meff_B ;
    Histo1DPtr _hist_meff_C ;
    Histo1DPtr _hist_meff_D ;
    Histo1DPtr _hist_meff_E ;

  };


  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1125961);

}
