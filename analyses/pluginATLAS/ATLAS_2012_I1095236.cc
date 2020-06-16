// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @author Peter Richardson
  class ATLAS_2012_I1095236 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1095236()
      : Analysis("ATLAS_2012_I1095236")
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
      declare(ChargedFinalState(Cuts::abseta < 3.0),"cfs");

      // Used for pTmiss
      declare(VisibleFinalState((Cuts::etaIn(-4.9,4.9))),"vfs");

      // Book histograms
      book(_count_SR0_A1 ,"count_SR0_A1", 1, 0., 1.);
      book(_count_SR0_B1 ,"count_SR0_B1", 1, 0., 1.);
      book(_count_SR0_C1 ,"count_SR0_C1", 1, 0., 1.);
      book(_count_SR0_A2 ,"count_SR0_A2", 1, 0., 1.);
      book(_count_SR0_B2 ,"count_SR0_B2", 1, 0., 1.);
      book(_count_SR0_C2 ,"count_SR0_C2", 1, 0., 1.);
      book(_count_SR1_D  ,"count_SR1_D" , 1, 0., 1.);
      book(_count_SR1_E  ,"count_SR1_E" , 1, 0., 1.);

      book(_hist_meff_SR0_A1   ,"hist_m_eff_SR0_A1", 14, 400., 1800.);
      book(_hist_meff_SR0_A2   ,"hist_m_eff_SR0_A2", 14, 400., 1800.);
      book(_hist_meff_SR1_D_e  ,"hist_meff_SR1_D_e" , 16, 600., 2200.);
      book(_hist_meff_SR1_D_mu ,"hist_meff_SR1_D_mu", 16, 600., 2200.);

      book(_hist_met_SR0_A1    ,"hist_met_SR0_A1", 14, 0., 700.);
      book(_hist_met_SR0_A2    ,"hist_met_SR0_A2", 14, 0., 700.);
      book(_hist_met_SR0_D_e   ,"hist_met_SR1_D_e" , 15, 0., 600.);
      book(_hist_met_SR0_D_mu  ,"hist_met_SR1_D_mu", 15, 0., 600.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      Jets cand_jets;
      const Jets jets = apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV);
      for (const Jet& jet : jets) {
        if ( fabs( jet.eta() ) < 2.8 ) {
          cand_jets.push_back(jet);
        }
      }

      const Particles cand_e  = apply<IdentifiedFinalState>(event, "elecs").particlesByPt();

      const Particles cand_mu = apply<IdentifiedFinalState>(event, "muons").particlesByPt();
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

      // get the loose leptons used to define the 0 lepton channel
      Particles loose_e, loose_mu;
      for ( const Particle & e  : cand_e ) {
        bool away = true;
        for ( const Jet& jet  : recon_jets ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away ) loose_e.push_back( e );
      }
      for ( const Particle & mu  : cand_mu ) {
        bool away = true;
        for ( const Jet& jet  : recon_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away ) loose_mu.push_back( mu );
      }
      // tight leptons for the 1-lepton channel
      Particles tight_mu;
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();
      for ( const Particle & mu  : loose_mu) {
        if(mu.perp()<20.) continue;
        double pTinCone = -mu.pT();
        for ( const Particle & track  : chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) <= 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 1.8*GeV )
          tight_mu.push_back(mu);
      }
      Particles tight_e;
      for ( const Particle & e  : loose_e ) {
        if(e.perp()<25.) continue;
        double pTinCone = -e.perp();
        for ( const Particle & track  : chg_tracks ) {
          if ( deltaR(e.momentum(),track.momentum()) <= 0.2 )
            pTinCone += track.pT();
        }
        if (pTinCone/e.perp()<0.1) {
          tight_e.push_back(e);
        }
      }

      // pTmiss
      Particles vfs_particles =
        apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle & p : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // get the number of b-tagged jets
      unsigned int ntagged=0;
      for (const Jet & jet : recon_jets ) {
        if(jet.perp()>50. && abs(jet.eta())<2.5 &&
           jet.bTagged() && rand()/static_cast<double>(RAND_MAX)<=0.60)
          ++ntagged;
      }

      // ATLAS calo problem
      if(rand()/static_cast<double>(RAND_MAX)<=0.42) {
        for ( const Jet & jet : recon_jets ) {
          double eta = jet.rapidity();
          double phi = jet.azimuthalAngle(MINUSPI_PLUSPI);
          if(jet.perp()>50 && eta>-0.1&&eta<1.5&&phi>-0.9&&phi<-0.5)
            vetoEvent;
        }
      }

      // at least 1 b tag
      if(ntagged==0) vetoEvent;

      // minumum Et miss
      if(eTmiss<80.) vetoEvent;

      // at least 3 jets pT > 50
      if(recon_jets.size()<3 || recon_jets[2].perp()<50.)
        vetoEvent;

      // m_eff
      double m_eff =  eTmiss;
      for(unsigned int ix=0;ix<3;++ix)
        m_eff += recon_jets[ix].perp();

      // delta Phi
      double min_dPhi = 999.999;
      double pTmiss_phi = pTmiss.phi();
      for(unsigned int ix=0;ix<3;++ix) {
        min_dPhi = min( min_dPhi, deltaPhi( pTmiss_phi, recon_jets[ix].phi() ) );
      }

      // 0-lepton channels
      if(loose_e.empty() && loose_mu.empty() &&
         recon_jets[0].perp()>130.  && eTmiss>130. &&
         eTmiss/m_eff>0.25 && min_dPhi>0.4) {
        // jet charge cut
        bool jetCharge = true;
        for(unsigned int ix=0;ix<3;++ix) {
          if(fabs(recon_jets[ix].eta())>2.) continue;
          double trackpT=0;
          for(const Particle & p : recon_jets[ix].particles()) {
            if(PID::charge3(p.pid())==0) continue;
            trackpT += p.perp();
          }
          if(trackpT/recon_jets[ix].perp()<0.05)
            jetCharge = false;
        }
        if(jetCharge) {
          // SR0-A region
          if(m_eff>500.) {
            _count_SR0_A1->fill(0.5,weight);
            _hist_meff_SR0_A1->fill(m_eff,weight);
            _hist_met_SR0_A1 ->fill(eTmiss,weight);
            if(ntagged>=2) {
              _count_SR0_A2->fill(0.5,weight);
              _hist_meff_SR0_A2->fill(m_eff,weight);
              _hist_met_SR0_A2 ->fill(eTmiss,weight);
            }
          }
          // SR0-B
          if(m_eff>700.) {
            _count_SR0_B1->fill(0.5,weight);
            if(ntagged>=2) _count_SR0_B2->fill(0.5,weight);
          }
          // SR0-C
          if(m_eff>900.) {
            _count_SR0_C1->fill(0.5,weight);
            if(ntagged>=2) _count_SR0_C2->fill(0.5,weight);
          }
        }
      }

      // 1-lepton channels
      if(tight_e.size() + tight_mu.size() == 1 &&
         recon_jets.size()>=4 && recon_jets[3].perp()>50.&&
         recon_jets[0].perp()>60.) {
        Particle lepton = tight_e.empty() ? tight_mu[0] : tight_e[0];
        m_eff += lepton.perp() + recon_jets[3].perp();
        // transverse mass cut
        double mT = 2.*(lepton.perp()*eTmiss-
                        lepton.px()*pTmiss.px()-
                        lepton.py()*pTmiss.py());
        mT = sqrt(mT);
        if(mT>100.&&m_eff>700.) {
          // D region
          _count_SR1_D->fill(0.5,weight);
          if(lepton.abspid()==PID::ELECTRON) {
            _hist_meff_SR1_D_e->fill(m_eff,weight);
            _hist_met_SR0_D_e->fill(eTmiss,weight);
          }
          else {
            _hist_meff_SR1_D_mu->fill(m_eff,weight);
            _hist_met_SR0_D_mu->fill(eTmiss,weight);
          }
          // E region
          if(eTmiss>200.) {
            _count_SR1_E->fill(0.5,weight);
          }
        }
      }
    }


    void finalize() {

      double norm = crossSection()/femtobarn*2.05/sumOfWeights();
      // these are number of events at 2.05fb^-1 per 100 GeV
      scale( _hist_meff_SR0_A1   , 100. * norm );
      scale( _hist_meff_SR0_A2   , 100. * norm );
      scale( _hist_meff_SR1_D_e  , 100. * norm );
      scale( _hist_meff_SR1_D_mu , 100. * norm );
      // these are number of events at 2.05fb^-1 per 50 GeV
      scale( _hist_met_SR0_A1, 50. * norm );
      scale( _hist_met_SR0_A2, 40. * norm );
      // these are number of events at 2.05fb^-1 per 40 GeV
      scale( _hist_met_SR0_D_e , 40. * norm );
      scale( _hist_met_SR0_D_mu, 40. * norm );
      // these are number of events at 2.05fb^-1
      scale(_count_SR0_A1,norm);
      scale(_count_SR0_B1,norm);
      scale(_count_SR0_C1,norm);
      scale(_count_SR0_A2,norm);
      scale(_count_SR0_B2,norm);
      scale(_count_SR0_C2,norm);
      scale(_count_SR1_D ,norm);
      scale(_count_SR1_E ,norm);
    }

    //@}


  private:

    Histo1DPtr _count_SR0_A1;
    Histo1DPtr _count_SR0_B1;
    Histo1DPtr _count_SR0_C1;
    Histo1DPtr _count_SR0_A2;
    Histo1DPtr _count_SR0_B2;
    Histo1DPtr _count_SR0_C2;
    Histo1DPtr _count_SR1_D;
    Histo1DPtr _count_SR1_E;

    Histo1DPtr _hist_meff_SR0_A1;
    Histo1DPtr _hist_meff_SR0_A2;
    Histo1DPtr _hist_meff_SR1_D_e;
    Histo1DPtr _hist_meff_SR1_D_mu;
    Histo1DPtr _hist_met_SR0_A1;
    Histo1DPtr _hist_met_SR0_A2;
    Histo1DPtr _hist_met_SR0_D_e;
    Histo1DPtr _hist_met_SR0_D_mu;

  };


  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1095236);

}
