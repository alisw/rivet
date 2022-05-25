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


  /// High-jet-multiplicity squark and gluino search
  class ATLAS_2011_S9225137 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2011_S9225137);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // veto region electrons
      Cut vetocut = Cuts::absetaIn(1.37, 1.52);
      IdentifiedFinalState veto_elecs(vetocut && Cuts::pT > 10*GeV);
      veto_elecs.acceptIdPair(PID::ELECTRON);
      declare(veto_elecs, "veto_elecs");

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 20*GeV);
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
      declare(ChargedFinalState(Cuts::abseta < 3), "cfs");

      /// Book histograms
      book(_etmisspT_55_NJ_6_obs , 1,1,1);
      book(_etmisspT_55_NJ_6_bac , 1,1,2);
      book(_etmisspT_55_NJ_6_sig , 1,1,3);
      book(_etmisspT_55_NJ_7_obs ,13,1,1);
      book(_etmisspT_55_NJ_7_bac ,13,1,2);
      book(_etmisspT_55_NJ_7_sig ,13,1,3);
      book(_etmisspT_55_NJ_8_obs ,15,1,1);
      book(_etmisspT_55_NJ_8_bac ,15,1,2);
      book(_etmisspT_55_NJ_8_sig ,15,1,3);
      book(_etmisspT_80_NJ_5_obs , 2,1,1);
      book(_etmisspT_80_NJ_5_bac , 2,1,2);
      book(_etmisspT_80_NJ_5_sig , 2,1,3);
      book(_etmisspT_80_NJ_6_obs ,14,1,1);
      book(_etmisspT_80_NJ_6_bac ,14,1,2);
      book(_etmisspT_80_NJ_6_sig ,14,1,3);
      book(_etmisspT_80_NJ_7_obs ,16,1,1);
      book(_etmisspT_80_NJ_7_bac ,16,1,2);
      book(_etmisspT_80_NJ_7_sig ,16,1,3);

      book(_njet55A_obs , 3,1,1);
      book(_njet55A_bac , 3,1,2);
      book(_njet55A_sig , 3,1,3);
      book(_njet55B_obs , 4,1,1);
      book(_njet55B_bac , 4,1,2);
      book(_njet55B_sig , 4,1,3);
      book(_njet55C_obs ,17,1,1);
      book(_njet55C_bac ,17,1,2);
      book(_njet55C_sig ,17,1,3);
      book(_njet80A_obs , 5,1,1);
      book(_njet80A_bac , 5,1,2);
      book(_njet80A_sig , 5,1,3);
      book(_njet80B_obs , 6,1,1);
      book(_njet80B_bac , 6,1,2);
      book(_njet80B_sig , 6,1,3);
      book(_njet80C_obs ,18,1,1);
      book(_njet80C_bac ,18,1,2);
      book(_njet80C_sig ,18,1,3);

      book(_count_7j55 ,"count_7j55", 1, 0., 1.);
      book(_count_8j55 ,"count_8j55", 1, 0., 1.);
      book(_count_6j80 ,"count_6j80", 1, 0., 1.);
      book(_count_7j80 ,"count_7j80", 1, 0., 1.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;
      // apply electron veto region
      Particles veto_e
        = apply<IdentifiedFinalState>(event, "veto_elecs").particles();
      if ( ! veto_e.empty() ) {
        MSG_DEBUG("electrons in veto region");
        vetoEvent;
      }

      // get the jet candidates
      Jets cand_jets;
      for (const Jet& jet :
               apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.eta() ) < 4.9 ) {
          cand_jets.push_back(jet);
        }
      }

      // candidate muons
      Particles cand_mu;
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();
      for ( const Particle & mu :
                apply<IdentifiedFinalState>(event, "muons").particlesByPt() ) {
        double pTinCone = -mu.pT();
        for ( const Particle & track : chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) <= 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 1.8*GeV )
          cand_mu.push_back(mu);
      }

      // candidate electrons

      Particles cand_e  =
        apply<IdentifiedFinalState>(event, "elecs").particlesByPt();

      // resolve jet/lepton ambiguity
      Jets cand_jets_2;
      for ( const Jet& jet : cand_jets ) {
        // candidates above eta=2.8 are jets
        if ( fabs( jet.eta() ) >= 2.8 )
          cand_jets_2.push_back( jet );
        // otherwise more the R=0.2 from an electrons
        else {
          bool away_from_e = true;
          for ( const Particle & e : cand_e ) {
            if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
              away_from_e = false;
              break;
            }
          }
          if ( away_from_e )
            cand_jets_2.push_back( jet );
        }
      }

      // only keep electrons more than R=0.4 from jets
      Particles recon_e;
      for ( const Particle & e : cand_e ) {
        bool away = true;
        for ( const Jet& jet : cand_jets_2 ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away )
          recon_e.push_back( e );
      }

      // only keep muons more than R=0.4 from jets
      Particles recon_mu;
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
      Particles vfs_particles =
        apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle & p : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // final jet filter
      Jets recon_jets;
      for (const Jet& jet : cand_jets_2) {
        if (jet.abseta() <= 2.8 ) recon_jets.push_back( jet );
      }

      // now only use recon_jets, recon_mu, recon_e

      // reject events with electrons and muons
      if ( ! ( recon_mu.empty() && recon_e.empty() ) ) {
        MSG_DEBUG("Charged leptons left after selection");
        vetoEvent;
      }

      // calculate H_T
      double HT = 0;
      for (const Jet& jet : recon_jets) {
        if (jet.pT() > 40*GeV) HT += jet.pT() ;
      }

      // number of jets and deltaR
      bool pass55DeltaR=true, pass80DeltaR=true;
      size_t njet55=0, njet80=0;
      for (size_t ix=0; ix < recon_jets.size(); ++ix) {
        if (recon_jets[ix].pT() > 80*GeV) ++njet80;
        if (recon_jets[ix].pT() > 55*GeV) ++njet55;
        for (size_t iy = ix + 1; iy < recon_jets.size(); ++iy) {
          if (recon_jets[ix].pT() > 55*GeV &&
              recon_jets[iy].pT() > 55*GeV &&
              deltaR(recon_jets[ix], recon_jets[iy]) < 0.6)
            pass55DeltaR = false;
          // if (recon_jets[ix].pT() > 80*GeV &&
          //     recon_jets[iy].pT() > 80*GeV &&
          //     deltaR(recon_jets[ix], recon_jets[iy]) < 0.6)
          //   pass80DeltaR = false;
        }
      }

      // require at least four jets with et > 55
      if (njet55 <= 3) vetoEvent;

      // plots of etmiss/ht
      double etht = eTmiss/sqrt(HT);
      if (njet55 == 6) {
        _etmisspT_55_NJ_6_obs->fill(etht,weight);
        _etmisspT_55_NJ_6_bac->fill(etht,weight);
        _etmisspT_55_NJ_6_sig->fill(etht,weight);
      } else if (njet55 == 7) {
        _etmisspT_55_NJ_7_obs->fill(etht,weight);
        _etmisspT_55_NJ_7_bac->fill(etht,weight);
        _etmisspT_55_NJ_7_sig->fill(etht,weight);
      } else if (njet55 == 8) {
        _etmisspT_55_NJ_8_obs->fill(etht,weight);
        _etmisspT_55_NJ_8_bac->fill(etht,weight);
        _etmisspT_55_NJ_8_sig->fill(etht,weight);
      }
      if (njet80 == 5) {
        _etmisspT_80_NJ_5_obs->fill(etht,weight);
        _etmisspT_80_NJ_5_bac->fill(etht,weight);
        _etmisspT_80_NJ_5_sig->fill(etht,weight);
      } else if (njet80 == 6) {
        _etmisspT_80_NJ_6_obs->fill(etht,weight);
        _etmisspT_80_NJ_6_bac->fill(etht,weight);
        _etmisspT_80_NJ_6_sig->fill(etht,weight);
      } else if (njet80 == 7) {
        _etmisspT_80_NJ_7_obs->fill(etht,weight);
        _etmisspT_80_NJ_7_bac->fill(etht,weight);
        _etmisspT_80_NJ_7_sig->fill(etht,weight);
      }

      if (etht > 1.5 && etht < 2.) {
        if (njet55 > 3) {
          _njet55A_obs->fill(njet55,weight);
          _njet55A_bac->fill(njet55,weight);
          _njet55A_sig->fill(njet55,weight);
        }
        if (njet80 > 3) {
          _njet80A_obs->fill(njet80,weight);
          _njet80A_bac->fill(njet80,weight);
          _njet80A_sig->fill(njet80,weight);
        }
      } else if (etht > 2. && etht < 3.) {
        if (njet55 > 3) {
          _njet55B_obs->fill(njet55,weight);
          _njet55B_bac->fill(njet55,weight);
          _njet55B_sig->fill(njet55,weight);
        }
        if (njet80 > 3) {
          _njet80B_obs->fill(njet80,weight);
          _njet80B_bac->fill(njet80,weight);
          _njet80B_sig->fill(njet80,weight);
        }
      } else {
        if (njet55 > 3) {
          _njet55C_obs->fill(njet55,weight);
          _njet55C_bac->fill(njet55,weight);
          _njet55C_sig->fill(njet55,weight);
        }
        if (njet80 > 3) {
          _njet80C_obs->fill(njet80,weight);
          _njet80C_bac->fill(njet80,weight);
          _njet80C_sig->fill(njet80,weight);
        }
      }

      // apply E_T/sqrt(H_T) cut
      if (etht <= 3.5*GeV) {
        MSG_DEBUG("Fails ET/sqrt(HT) cut ");
        vetoEvent;
      }

      // check passes at least one delta5/ njet number cut
      if (!(pass55DeltaR && njet55 >= 7) && !(pass80DeltaR && njet80 >= 6) ) {
        MSG_DEBUG("Fails DeltaR cut or jet number cuts");
        vetoEvent;
      }

      // 7j55
      if (njet55 >= 7 && pass55DeltaR) _count_7j55->fill( 0.5, weight);
      // 8j55
      if (njet55 >= 8 && pass55DeltaR) _count_8j55->fill( 0.5, weight);
      // 6j80
      if (njet80 >= 6 && pass80DeltaR) _count_6j80->fill( 0.5, weight);
      // 7j80
      if (njet80 >= 7 && pass80DeltaR) _count_7j80->fill( 0.5, weight);
    }

    //@}

    void finalize() {
      const double norm = crossSection()/femtobarn*1.34/sumOfWeights();
      scale(_etmisspT_55_NJ_6_obs,norm);
      scale(_etmisspT_55_NJ_6_bac,norm);
      scale(_etmisspT_55_NJ_6_sig,norm);
      scale(_etmisspT_55_NJ_7_obs,norm);
      scale(_etmisspT_55_NJ_7_bac,norm);
      scale(_etmisspT_55_NJ_7_sig,norm);
      scale(_etmisspT_55_NJ_8_obs,norm);
      scale(_etmisspT_55_NJ_8_bac,norm);
      scale(_etmisspT_55_NJ_8_sig,norm);
      scale(_etmisspT_80_NJ_5_obs,norm);
      scale(_etmisspT_80_NJ_5_bac,norm);
      scale(_etmisspT_80_NJ_5_sig,norm);
      scale(_etmisspT_80_NJ_6_obs,norm);
      scale(_etmisspT_80_NJ_6_bac,norm);
      scale(_etmisspT_80_NJ_6_sig,norm);
      scale(_etmisspT_80_NJ_7_obs,norm);
      scale(_etmisspT_80_NJ_7_bac,norm);
      scale(_etmisspT_80_NJ_7_sig,norm);
      scale(_njet55A_obs,norm);
      scale(_njet55A_bac,norm);
      scale(_njet55A_sig,norm);
      scale(_njet55B_obs,norm);
      scale(_njet55B_bac,norm);
      scale(_njet55B_sig,norm);
      scale(_njet55C_obs,norm);
      scale(_njet55C_bac,norm);
      scale(_njet55C_sig,norm);
      scale(_njet80A_obs,norm);
      scale(_njet80A_bac,norm);
      scale(_njet80A_sig,norm);
      scale(_njet80B_obs,norm);
      scale(_njet80B_bac,norm);
      scale(_njet80B_sig,norm);
      scale(_njet80C_obs,norm);
      scale(_njet80C_bac,norm);
      scale(_njet80C_sig,norm);
      scale(_count_7j55,norm);
      scale(_count_8j55,norm);
      scale(_count_6j80,norm);
      scale(_count_7j80,norm);
    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _etmisspT_55_NJ_6_obs;
    Histo1DPtr _etmisspT_55_NJ_6_bac;
    Histo1DPtr _etmisspT_55_NJ_6_sig;
    Histo1DPtr _etmisspT_55_NJ_7_obs;
    Histo1DPtr _etmisspT_55_NJ_7_bac;
    Histo1DPtr _etmisspT_55_NJ_7_sig;
    Histo1DPtr _etmisspT_55_NJ_8_obs;
    Histo1DPtr _etmisspT_55_NJ_8_bac;
    Histo1DPtr _etmisspT_55_NJ_8_sig;
    Histo1DPtr _etmisspT_80_NJ_5_obs;
    Histo1DPtr _etmisspT_80_NJ_5_bac;
    Histo1DPtr _etmisspT_80_NJ_5_sig;
    Histo1DPtr _etmisspT_80_NJ_6_obs;
    Histo1DPtr _etmisspT_80_NJ_6_bac;
    Histo1DPtr _etmisspT_80_NJ_6_sig;
    Histo1DPtr _etmisspT_80_NJ_7_obs;
    Histo1DPtr _etmisspT_80_NJ_7_bac;
    Histo1DPtr _etmisspT_80_NJ_7_sig;
    Histo1DPtr _njet55A_obs;
    Histo1DPtr _njet55A_bac;
    Histo1DPtr _njet55A_sig;
    Histo1DPtr _njet55B_obs;
    Histo1DPtr _njet55B_bac;
    Histo1DPtr _njet55B_sig;
    Histo1DPtr _njet55C_obs;
    Histo1DPtr _njet55C_bac;
    Histo1DPtr _njet55C_sig;
    Histo1DPtr _njet80A_obs;
    Histo1DPtr _njet80A_bac;
    Histo1DPtr _njet80A_sig;
    Histo1DPtr _njet80B_obs;
    Histo1DPtr _njet80B_bac;
    Histo1DPtr _njet80B_sig;
    Histo1DPtr _njet80C_obs;
    Histo1DPtr _njet80C_bac;
    Histo1DPtr _njet80C_sig;
    Histo1DPtr _count_7j55;
    Histo1DPtr _count_8j55;
    Histo1DPtr _count_6j80;
    Histo1DPtr _count_7j80;
    //@}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ATLAS_2011_S9225137, ATLAS_2011_I939504);

}
