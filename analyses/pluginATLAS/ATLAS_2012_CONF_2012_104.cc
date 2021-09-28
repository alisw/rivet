// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2012_CONF_2012_104 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_CONF_2012_104()
      : Analysis("ATLAS_2012_CONF_2012_104")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialize projections before the run
    void init() {

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 10*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // Jet finder
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3 && Cuts::pT > 0.5*GeV), "cfs");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9),"vfs");

      // Book histograms
      book(_count_e  ,"count_e" , 1, 0., 1.);
      book(_count_mu ,"count_mu", 1, 0., 1.);

      book(_hist_eTmiss_e  ,"hist_eTmiss_e"  , 25, 0., 1000.);
      book(_hist_eTmiss_mu ,"hist_eTmiss_mu" , 25, 0., 1000.);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // get the candiate jets
      Jets cand_jets;
      for ( const Jet& jet :
                apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.eta() ) < 2.8 ) {
          cand_jets.push_back(jet);
        }
      }

      // get the candidate "medium" leptons without isolation
      Particles cand_e;
      for( const Particle & e :
               apply<IdentifiedFinalState>(event, "elecs").particlesByPt()) {
        // remove any leptons within 0.4 of any candidate jets
        bool e_near_jet = false;
        for ( const Jet& jet : cand_jets ) {
          double dR = deltaR(e.momentum(),jet.momentum());
          if ( dR < 0.4 && dR > 0.2 ) {
            e_near_jet = true;
            break;
          }
        }
        if ( ! e_near_jet ) cand_e.push_back(e);
      }
      Particles cand_mu;
      for( const Particle & mu :
               apply<IdentifiedFinalState>(event, "muons").particlesByPt()) {
        // remove any leptons within 0.4 of any candidate jets
        bool mu_near_jet = false;
        for ( const Jet& jet : cand_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            mu_near_jet = true;
            break;
          }
        }
        if ( ! mu_near_jet ) cand_mu.push_back(mu);
      }
      // apply the isolation
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();
      // pTcone around muon track (hard)
      Particles recon_mu;
      for ( const Particle & mu : cand_mu ) {
        double pTinCone = -mu.pT();
        if(-pTinCone<25.) continue;
        for ( const Particle & track : chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 1.8*GeV ) recon_mu.push_back(mu);
      }
      // pTcone around electron track (hard)
      Particles recon_e;
      for ( const Particle & e : cand_e ) {
        double pTinCone = -e.pT();
        if(-pTinCone<25.) continue;
        for ( const Particle & track : chg_tracks ) {
          if ( deltaR(e.momentum(),track.momentum()) < 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 0.1 * e.pT() ) recon_e.push_back(e);
      }

      // discard jets that overlap with electrons
      Jets recon_jets;
      for ( const Jet& jet : cand_jets ) {
        if(jet.abseta()>2.5||
           jet.perp()<25.) continue;
        bool away_from_e = true;
        for ( const Particle & e : cand_e ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e ) recon_jets.push_back( jet );
      }

      // pTmiss
      FourMomentum pTmiss;
      for ( const Particle & p :
                apply<VisibleFinalState>(event, "vfs").particles() ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // at least 4 jets with pT>80.
      if(recon_jets.size()<4 || recon_jets[3].perp()<80.) vetoEvent;

      // only 1 signal lepton
      if( recon_e.size() + recon_mu.size() != 1 )
        vetoEvent;
      if( cand_e .size() + cand_mu .size() != 1 )
        vetoEvent;

      // start of meff calculation
      double HT=0.;
      for( const Jet & jet : recon_jets) {
        double pT = jet.perp();
        if(pT>40.) HT += pT;
      }

      // get the lepton
      Particle lepton = recon_e.empty() ? recon_mu[0] : recon_e[0];

      // lepton variables
      double pT = lepton.perp();

      double mT  = 2.*(pT*eTmiss -
                       lepton.px()*pTmiss.px() -
                       lepton.py()*pTmiss.py());
      mT = sqrt(mT);
      HT += pT;
      double m_eff_inc  = HT + eTmiss + pT;
      double m_eff_4 = eTmiss + pT;
      for(unsigned int ix=0;ix<4;++ix)
        m_eff_4 +=  recon_jets[ix].perp();

      // four jet selecton
      if(mT>100.&& eTmiss/m_eff_4>0.2 &&
         m_eff_inc > 800.) {
        if( eTmiss > 250. ) {
          if(lepton.abspid()==PID::ELECTRON)
            _count_e->fill(0.5,weight);
          else if(lepton.abspid()==PID::MUON)
            _count_mu->fill(0.5,weight);
        }
        if(lepton.abspid()==PID::ELECTRON)
          _hist_eTmiss_e ->fill(eTmiss,weight);
        else if(lepton.abspid()==PID::MUON)
          _hist_eTmiss_mu->fill(eTmiss,weight);
      }
    }
    //@}


    void finalize() {

      double norm = 5.8* crossSection()/sumOfWeights()/femtobarn;
      scale(_count_e ,norm);
      scale(_count_mu,norm);
      scale(_hist_eTmiss_e  ,40.*norm);
      scale(_hist_eTmiss_mu ,40.*norm);

    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _count_e ;
    Histo1DPtr _count_mu;

    Histo1DPtr _hist_eTmiss_e ;
    Histo1DPtr _hist_eTmiss_mu;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_CONF_2012_104);

}
