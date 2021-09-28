// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  /// @author Chris Wymant
  class ATLAS_2011_S9212183 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_S9212183()
      : Analysis("ATLAS_2011_S9212183")
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
      declare(FastJets(FinalState(), FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // All tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3.0), "cfs");

      // Used for pTmiss (N.B. the real 'vfs' extends beyond 4.5 to |eta| = 4.9)
      declare(VisibleFinalState(Cuts::abseta < 4.5), "vfs");


      // Book histograms
      book(_count_2j   ,"count_2j",   1, 0., 1.);
      book(_count_3j   ,"count_3j",   1, 0., 1.);
      book(_count_4j5  ,"count_4j5",  1, 0., 1.);
      book(_count_4j10 ,"count_4j10", 1, 0., 1.);
      book(_count_HM   ,"count_HM",   1, 0., 1.);

      book(_hist_meff_2j  ,1, 1, 1);
      book(_hist_meff_3j  ,2, 1, 1);
      book(_hist_meff_4j  ,3, 1, 1);
      book(_hist_meff_HM  ,4, 1, 1);

      book(_hist_eTmiss  ,"Et_miss", 20, 0., 1000.);
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

      // Muon isolation not mentioned in hep-exp 1109.6572, unlike in 1102.5290,
      // but assumed to still be applicable
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
      Jets cand_jets_2;
      for ( const Jet& jet : cand_jets ) {
        if ( fabs( jet.eta() ) >= 2.8 ) {
          cand_jets_2.push_back( jet );
        } else {
          bool away_from_e = true;
          for ( const Particle & e : cand_e ) {
            if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
              away_from_e = false;
              break;
            }
          }
          if ( away_from_e ) cand_jets_2.push_back( jet );
        }
      }


      Particles recon_e, recon_mu;

      for ( const Particle & e : cand_e ) {
        bool away = true;
        for ( const Jet& jet : cand_jets_2 ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away ) recon_e.push_back( e );
      }

      for ( const Particle & mu : cand_mu ) {
        bool away = true;
        for ( const Jet& jet : cand_jets_2 ) {
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
      for ( const Jet& jet : cand_jets_2 ) {
        if ( fabs( jet.eta() ) > 4.5 ) pTmiss -= jet.momentum();
      }
      double eTmiss = pTmiss.pT();


      // Final jet filter
      Jets recon_jets;
      for ( const Jet& jet : cand_jets_2 ) {
        if ( fabs( jet.eta() ) <= 2.8 ) recon_jets.push_back( jet );
      }
      // NB. It seems that jets with |eta| > 2.8 could have been thrown away at
      // the start; we don't do so, in order to follow both the structure of
      // the paper and the similar Rivet analysis ATLAS_2011_S8983313

      // 'candidate' muons needed only 10 GeV, to cause a veto they need 20 GeV
      Particles veto_mu;
      for ( const Particle & mu : cand_mu ) {
        if ( mu.pT() >= 20.0*GeV ) veto_mu.push_back(mu);
      }

      if ( ! ( veto_mu.empty() && recon_e.empty() ) ) {
        MSG_DEBUG("Charged leptons left after selection");
        vetoEvent;
      }

      if ( eTmiss <= 130 * GeV ) {
        MSG_DEBUG("Not enough eTmiss: " << eTmiss << " < 130");
        vetoEvent;
      }

      if ( recon_jets.empty() || recon_jets[0].pT() <= 130.0 * GeV ) {
        MSG_DEBUG("No hard leading jet in " << recon_jets.size() << " jets");
        vetoEvent;
      }

      // ==================== observables ====================

      int Njets = 0;
      double min_dPhi = 999.999;
      double pTmiss_phi = pTmiss.phi();
      for ( const Jet& jet : recon_jets ) {
        if ( jet.pT() > 40 * GeV ) {
          if ( Njets < 3 ) {
            min_dPhi = min( min_dPhi, deltaPhi( pTmiss_phi, jet.phi() ) );
          }
          ++Njets;
        }
      }

      int NjetsHighMass = 0;
      for ( const Jet& jet : recon_jets ) {
        if ( jet.pT() > 80.0 * GeV ) {
          ++NjetsHighMass;
        }
      }

      if ( Njets < 2 ) {
        MSG_DEBUG("Only " << Njets << " >40 GeV jets left");
        vetoEvent;
      }

      if ( min_dPhi <= 0.4 ) {
        MSG_DEBUG("dPhi too small");
        vetoEvent;
      }

      // m_eff
      double m_eff_2j = eTmiss + recon_jets[0].pT() + recon_jets[1].pT();
      double m_eff_3j = recon_jets.size() < 3 ? -999.0 : m_eff_2j + recon_jets[2].pT();
      double m_eff_4j = recon_jets.size() < 4 ? -999.0 : m_eff_3j + recon_jets[3].pT();
      double m_eff_HM = eTmiss;
      for ( const Jet& jet : recon_jets ) {
        if ( jet.pT() > 40.0 * GeV ) m_eff_HM += jet.pT();
      }

      double et_meff_2j = eTmiss / m_eff_2j;
      double et_meff_3j = eTmiss / m_eff_3j;
      double et_meff_4j = eTmiss / m_eff_4j;
      double et_meff_HM = eTmiss / m_eff_HM;


      // ==================== FILL ====================

      MSG_DEBUG( "Trying to fill "
                 << Njets << ' '
                 << m_eff_2j << ' '
                 << et_meff_2j << ' '
                 << m_eff_3j << ' '
                 << et_meff_3j << ' '
                 << m_eff_4j << ' '
                 << et_meff_4j << ' '
                 << m_eff_HM << ' '
                 << et_meff_HM );


      _hist_eTmiss->fill(eTmiss, weight);


      // 2j region
      if ( et_meff_2j > 0.3 ) {
        _hist_meff_2j->fill(m_eff_2j, weight);
        if ( m_eff_2j > 1000 * GeV ) {
          MSG_DEBUG("Hits 2j");
          _count_2j->fill(0.5, weight);
        }
      }


      // 3j region
      if ( Njets >= 3 && et_meff_3j > 0.25 ) {
        _hist_meff_3j->fill(m_eff_3j, weight);
        if ( m_eff_3j > 1000 * GeV ) {
          MSG_DEBUG("Hits 3j");
          _count_3j->fill(0.5, weight);
        }
      }


      // 4j5 & 4j10 regions
      if ( Njets >= 4 && et_meff_4j > 0.25 ) {
        _hist_meff_4j->fill(m_eff_4j, weight);
        if ( m_eff_4j > 500 * GeV ) {
          MSG_DEBUG("Hits 4j5");
          _count_4j5->fill(0.5, weight);
        }
        if ( m_eff_4j > 1000 * GeV ) {
          MSG_DEBUG("Hits 4j10");
          _count_4j10->fill(0.5, weight);
        }
      }


      // High mass region
      if ( NjetsHighMass >= 4 && et_meff_HM > 0.2 ) {
        _hist_meff_HM->fill(m_eff_HM, weight);
        if ( m_eff_HM > 1100 * GeV ) {
          MSG_DEBUG("Hits HM");
          _count_HM->fill(0.5, weight);
        }
      }

    }


    void finalize() {
      // Two, three and four jet channels have bin width = 100 (GeV)
      // High mass channel has bin width = 150 (GeV)
      // Integrated luminosity = 1040 (pb)
      scale( _hist_meff_2j, 100. * 1040 * crossSection()/sumOfWeights() );
      scale( _hist_meff_3j, 100. * 1040 * crossSection()/sumOfWeights() );
      scale( _hist_meff_4j, 100. * 1040 * crossSection()/sumOfWeights() );
      scale( _hist_meff_HM, 150. * 1040 * crossSection()/sumOfWeights() );
    }

    //@}


  private:

    Histo1DPtr _count_2j;
    Histo1DPtr _count_3j;
    Histo1DPtr _count_4j5;
    Histo1DPtr _count_4j10;
    Histo1DPtr _count_HM;
    Histo1DPtr _hist_meff_2j;
    Histo1DPtr _hist_meff_3j;
    Histo1DPtr _hist_meff_4j;
    Histo1DPtr _hist_meff_HM;
    Histo1DPtr _hist_eTmiss;

  };


  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_S9212183);

}
