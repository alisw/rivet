// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetMT2.hh"

namespace Rivet {


  class ATLAS_2012_CONF_2012_153 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2012_CONF_2012_153()
      : Analysis("ATLAS_2012_CONF_2012_153")
    {    }

    //@}



  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 10*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");


      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9), "vfs");

      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);

      /// Jet finder
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3.0), "cfs");

      vector<double> edges_meff;
      edges_meff.push_back(   0);
      edges_meff.push_back( 150);
      edges_meff.push_back( 300);
      edges_meff.push_back( 500);
      edges_meff.push_back(1000);
      edges_meff.push_back(1500);

      vector<double> edges_eT;
      edges_eT.push_back(0);
      edges_eT.push_back(50);
      edges_eT.push_back(150);
      edges_eT.push_back(300);
      edges_eT.push_back(500);

      // Book histograms
      book(_hist_electrons ,"hist_electrons_before", 11, -0.5,10.5);
      book(_hist_muons     ,"hist_muons_before"    , 11, -0.5,10.5);
      book(_hist_leptons   ,"hist_leptons_before"  , 11, -0.5,10.5);
      book(_hist_4leptons  ,"hist_4leptons", 1, 0.,1.);
      book(_hist_veto      ,"hist_veto", 1, 0., 1.);
      book(_hist_etmiss    ,"hist_etmiss",edges_eT);
      book(_hist_meff      ,"hist_m_eff",edges_meff);
      book(_count_SR1      ,"count_SR1", 1, 0., 1.);
      book(_count_SR2      ,"count_SR2", 1, 0., 1.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;
      // get the jet candidates
      Jets cand_jets;
      for (const Jet& jet : apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if (jet.abseta() < 2.5) cand_jets.push_back(jet);
      }

      // candidate muons
      Particles cand_mu = apply<IdentifiedFinalState>(event, "muons").particlesByPt();

      // candidate electrons
      // Discard if two electrons are within R=0.1
      Particles temp = apply<IdentifiedFinalState>(event, "elecs").particles(cmpMomByE);
      vector<bool> vetoed(temp.size(),false);
      Particles cand_e;
      for (size_t ix = 0; ix < temp.size(); ++ix) {
        if (vetoed[ix]) continue;
        for (size_t iy = ix+1; iy < temp.size(); ++iy) {
          if ( deltaR(temp[ix], temp[iy]) < 0.1 ) vetoed[iy] = true;
        }
        if (!vetoed[ix]) cand_e.push_back(temp[ix]);
      }

      // Sort by transverse momentum
      sortByPt(cand_e);

      // resolve jet/lepton ambiguity
      Jets recon_jets;
      for ( const Jet& jet : cand_jets ) {
        bool away_from_e = true;
        for ( const Particle& e : cand_e ) {
          if (deltaR(e, jet) <= 0.2) {
            away_from_e = false;
            break;
          }
        }
        if (away_from_e) recon_jets.push_back( jet );
      }

      // only keep electrons more than R=0.4 from jets
      Particles cand2_e;
      for (const Particle& e : cand_e) {
        // at least 0.4 from any jets
        bool away = true;
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(e, jet) < 0.4 ) {
            away = false;
            break;
          }
        }
        // if isolated keep it
        if ( away )
          cand2_e.push_back( e );
      }

      // only keep muons more than R=0.4 from jets
      Particles cand2_mu;
      for(const Particle & mu : cand_mu ) {
        bool away = true;
        // at least 0.4 from any jets
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(mu, jet) < 0.4 ) {
            away = false;
            break;
          }
        }
        if (away) cand2_mu.push_back( mu );
      }

      // electron and muon more than 0.1 apart
      Particles cand3_e;
      for ( const Particle & e : cand2_e ) {
        bool away = true;
        for( const Particle & mu : cand2_mu ) {
          if( deltaR(e, mu) < 0.1) {
            away = false;
            break;
          }
        }
        if (away) cand3_e.push_back(e);
      }
      Particles cand3_mu;
      for( const Particle & mu : cand2_mu ) {
        bool away = true;
        for ( const Particle & e : cand2_e ) {
          if( deltaR(e, mu) < 0.1) {
            away = false;
            break;
          }
        }
        if (away) cand3_mu.push_back(mu);
      }

      // pTmiss
      Particles vfs_particles =
        apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle & p : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // apply electron isolation
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();
      Particles cand4_e;
      for (const Particle& e : cand3_e) {
        // charge isolation
        double pTinCone = -e.pT();
        for (const Particle& track : chg_tracks) {
          if (track.pT() > 0.4*GeV && deltaR(e, track) <= 0.3 )
            pTinCone += track.pT();
        }
        if (pTinCone/e.pT() > 0.16) continue;
        // all particles isolation
        pTinCone = -e.pT();
        for (const Particle& p : vfs_particles) {
          if (p.abspid() != PID::MUON && deltaR(e, p) <= 0.3 )
            pTinCone += p.pT();
        }
        if (pTinCone/e.pT() < 0.18) cand4_e.push_back(e);
      }

      // apply muon isolation
      Particles cand4_mu;
      for ( const Particle & mu : cand3_mu ) {
        double pTinCone = -mu.perp();
        for ( const Particle & track : chg_tracks ) {
          if (track.pT() > 1*GeV && deltaR(mu, track) <= 0.3)
            pTinCone += track.pT();
        }
        if (pTinCone/mu.pT() < 0.12) cand4_mu.push_back(mu);
      }

      // same SOSF pairs m>12.
      Particles recon_e;
      for(const Particle& e : cand4_e) {
        bool veto = false;
        for(const Particle& e2 : cand4_e) {
          if (e.pid()*e2.pid() < 0 && (e.momentum()+e2.momentum()).mass() < 12*GeV) {
            veto = true;
            break;
          }
        }
        if (!veto) recon_e.push_back(e);
      }
      Particles recon_mu;
      for(const Particle& mu : cand4_mu) {
        bool veto = false;
        for(const Particle& mu2 : cand4_mu) {
          if (mu.pid()*mu2.pid() < 0 && (mu.momentum()+mu2.momentum()).mass() < 12*GeV) {
            veto = true;
            break;
          }
        }
        if (!veto) recon_mu.push_back(mu);
      }

      // now only use recon_jets, recon_mu, recon_e
      _hist_electrons->fill(recon_e.size(), weight);
      _hist_muons->fill(recon_mu.size(), weight);
      _hist_leptons->fill(recon_mu.size() + recon_e.size(), weight);
      if (recon_mu.size() + recon_e.size() > 3) {
        _hist_4leptons->fill(0.5, weight);
      }

      // reject events with less than 4 electrons and muons
      if (recon_mu.size() + recon_e.size() < 4) {
        MSG_DEBUG("To few charged leptons left after selection");
        vetoEvent;
      }


      // or two lepton trigger
      bool passDouble =
        (recon_mu.size()>=2 && ( (recon_mu[1].pT()>14*GeV) ||
                                 (recon_mu[0].pT()>18*GeV && recon_mu[1].perp() > 10*GeV) )) ||
        (recon_e.size() >=2 && ( (recon_e [1].pT()>14*GeV) ||
                                 (recon_e [0].pT()>25*GeV && recon_e [1].perp() > 10*GeV) )) ||
        (!recon_e.empty() && !recon_mu.empty() &&
         ( (recon_e[0].pT() > 14*GeV && recon_mu[0].pT() > 10*GeV)||
           (recon_e[0].pT() > 10*GeV && recon_mu[0].pT() > 18*GeV) ));

      // must pass a trigger
      if (!passDouble ) {
        MSG_DEBUG("Hardest lepton fails trigger");
        _hist_veto->fill(0.5, weight);
        vetoEvent;
      }

      // calculate meff
      double meff = eTmiss;
      for ( const Particle & e  : recon_e  ) meff += e.perp();
      for ( const Particle & mu : recon_mu ) meff += mu.perp();
      for ( const Jet & jet : recon_jets ) {
        const double pT = jet.pT();
        if (pT > 40*GeV) meff += pT;
      }

      // 2/3 leptons --> find 1 SFOS pair in range and veto event
      // 4+  leptons --> find 2 SFOS pairs and in range veto event
      for (size_t ix = 0; ix < recon_e.size(); ++ix) {
        for (size_t iy = ix+1; iy < recon_e.size(); ++iy) {
          if (recon_e[ix].pid()*recon_e[iy].pid() > 0) continue;
          const FourMomentum ppair = recon_e[ix].momentum() + recon_e[iy].momentum();
          if (inRange(ppair.mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;

          // check triplets with electron
          for (size_t iz = 0; iz < recon_e.size(); ++iz) {
            if (iz == ix || iz == iy) continue;
            if (inRange((ppair+recon_e[iz].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
          }

          // check triplets with muon
          for (size_t iz = 0; iz < recon_mu.size(); ++iz) {
            if (inRange((ppair+recon_mu[iz].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
          }

          // check quadruplets with electrons
          for (size_t iz = 0; iz < recon_e.size(); ++iz) {
            for (size_t iw = iz+1; iw < recon_e.size(); ++iw) {
              if (iz==ix || iz==iy || iw==ix || iw==iy) continue;
              if (recon_e[iz].pid()*recon_e[iw].pid() > 0) continue;
              if (inRange((ppair+recon_e[iz].momentum()+recon_e[iw].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
            }
          }
          // check quadruplets with muons
          for (size_t iz = 0; iz < recon_mu.size(); ++iz) {
            for (size_t iw = iz+1; iw < recon_mu.size(); ++iw) {
              if (recon_mu[iz].pid()*recon_mu[iw].pid() > 0) continue;
              if (inRange((ppair+recon_mu[iz].momentum()+recon_mu[iw].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
            }
          }
        }
      }

      // Muon pairs
      for (size_t ix = 0; ix < recon_mu.size(); ++ix) {
        for (size_t iy = ix+1; iy < recon_mu.size(); ++iy) {
          if (recon_mu[ix].pid()*recon_mu[iy].pid()>0) continue;
          const FourMomentum ppair = recon_mu[ix].momentum()+recon_mu[iy].momentum();
          if (inRange(ppair.mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;

          // check triplets with muon
          for (size_t iz = 0; iz < recon_mu.size(); ++iz) {
            if (iz==ix || iz==iy) continue;
            if (inRange((ppair+recon_mu[iz].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
          }

          // check triplets with electron
          for (size_t iz = 0; iz < recon_e.size(); ++iz) {
            if (inRange((ppair+recon_e[iz].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
          }

          // check muon quadruplets
          for (size_t iz = 0; iz < recon_mu.size(); ++iz) {
            for (size_t iw = iz+1; iy < recon_mu.size(); ++iy) {
              if (iz==ix || iz==iy || iw==ix || iw==iy) continue;
              if (recon_mu[iz].pid()*recon_mu[iw].pid() > 0) continue;
              if (inRange((ppair+recon_mu[iz].momentum()+recon_mu[iw].momentum()).mass(), 81.2*GeV, 101.2*GeV)) vetoEvent;
            }
          }
        }
      }

      // Make the control plots
      _hist_etmiss->fill(eTmiss,weight);
      _hist_meff  ->fill(meff  ,weight);
      // Finally the counts
      if (eTmiss > 50*GeV) _count_SR1->fill(0.5,weight);
      if (meff  >0*GeV) _count_SR2->fill(0.5,weight);

    }

    //@}

    void finalize() {
      double norm = crossSection()/femtobarn*13./sumOfWeights();
      scale(_hist_etmiss,norm*20.);
      scale(_hist_meff  ,norm*20.);
      scale(_count_SR1,norm);
      scale(_count_SR2,norm);
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hist_electrons;
    Histo1DPtr _hist_muons;
    Histo1DPtr _hist_leptons;
    Histo1DPtr _hist_4leptons;
    Histo1DPtr _hist_veto;
    Histo1DPtr _hist_etmiss;
    Histo1DPtr _hist_meff;
    Histo1DPtr _count_SR1;
    Histo1DPtr _count_SR2;
    //@}

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2012_CONF_2012_153);

}
