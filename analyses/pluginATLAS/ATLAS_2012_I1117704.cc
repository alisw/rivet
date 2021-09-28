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


  class ATLAS_2012_I1117704 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2012_I1117704()
      : Analysis("ATLAS_2012_I1117704")
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

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9), "vfs");

      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);

      /// Jet finder
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3),"cfs");

      /// Book histograms
      book(_etmiss_HT_7j55 ,"etmiss_HT_7j55", 8, 0., 16.);
      book(_etmiss_HT_8j55 ,"etmiss_HT_8j55", 8, 0., 16.);
      book(_etmiss_HT_9j55 ,"etmiss_HT_9j55", 8, 0., 16.);
      book(_etmiss_HT_6j80 ,"etmiss_HT_6j80", 8, 0., 16.);
      book(_etmiss_HT_7j80 ,"etmiss_HT_7j80", 8, 0., 16.);
      book(_etmiss_HT_8j80 ,"etmiss_HT_8j80", 8, 0., 16.);

      book(_hist_njet55 ,"hist_njet55", 11, 2.5, 13.5);
      book(_hist_njet80 ,"hist_njet80", 11, 2.5, 13.5);

      book(_count_7j55 ,"count_7j55", 1, 0., 1.);
      book(_count_8j55 ,"count_8j55", 1, 0., 1.);
      book(_count_9j55 ,"count_9j55", 1, 0., 1.);
      book(_count_6j80 ,"count_6j80", 1, 0., 1.);
      book(_count_7j80 ,"count_7j80", 1, 0., 1.);
      book(_count_8j80 ,"count_8j80", 1, 0., 1.);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      // get the jet candidates
      Jets cand_jets;
      for (const Jet& jet :
               apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.eta() ) < 2.8 ) {
          cand_jets.push_back(jet);
        }
      }

      // candidate muons
      Particles cand_mu;
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();
      for ( const Particle& mu :
                apply<IdentifiedFinalState>(event, "muons").particlesByPt() ) {
        double pTinCone = -mu.pT();
        for ( const Particle& track : chg_tracks ) {
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
      Jets recon_jets;
      for ( const Jet& jet : cand_jets ) {
        // candidates after |eta| < 2.8
        if ( fabs( jet.eta() ) >= 2.8 ) continue;
        bool away_from_e = true;
        for ( const Particle& e : cand_e ) {
          if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e ) recon_jets.push_back( jet );
      }

      // only keep electrons more than R=0.4 from jets
      Particles recon_e;
      for ( const Particle& e : cand_e ) {
        bool away = true;
        for ( const Jet& jet : recon_jets ) {
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
      for ( const Particle& mu : cand_mu ) {
        bool away = true;
        for ( const Jet& jet : recon_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if ( away )
          recon_mu.push_back( mu );
      }

      // pTmiss
      Particles vfs_particles =
        apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      for ( const Particle& p  : vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // now only use recon_jets, recon_mu, recon_e

      // reject events with electrons and muons
      if ( ! ( recon_mu.empty() && recon_e.empty() ) ) {
        MSG_DEBUG("Charged leptons left after selection");
        vetoEvent;
      }

      // calculate H_T
      double HT=0;
      for ( const Jet& jet : recon_jets ) {
        if ( jet.pT() > 40 * GeV )
          HT += jet.pT() ;
      }

      // number of jets
      unsigned int njet55=0, njet80=0;
      for (unsigned int ix=0;ix<recon_jets.size();++ix) {
        if(recon_jets[ix].pT()>80.*GeV) ++njet80;
        if(recon_jets[ix].pT()>55.*GeV) ++njet55;
      }

      if(njet55==0) vetoEvent;

      double ratio = eTmiss/sqrt(HT);

      if(ratio>4.) {
        _hist_njet55->fill(njet55,weight);
        _hist_njet80->fill(njet80,weight);
        // 7j55
        if(njet55>=7)
          _count_7j55->fill( 0.5, weight);
        // 8j55
        if(njet55>=8)
          _count_8j55->fill( 0.5, weight) ;
        // 8j55
        if(njet55>=9)
          _count_9j55->fill( 0.5, weight) ;
        // 6j80
        if(njet80>=6)
          _count_6j80->fill( 0.5, weight) ;
        // 7j80
        if(njet80>=7)
          _count_7j80->fill( 0.5, weight) ;
        // 8j80
        if(njet80>=8)
          _count_8j80->fill( 0.5, weight) ;
      }

      if(njet55>=7)
        _etmiss_HT_7j55->fill( ratio, weight);
      // 8j55
      if(njet55>=8)
        _etmiss_HT_8j55->fill( ratio, weight) ;
      // 8j55
      if(njet55>=9)
        _etmiss_HT_9j55->fill( ratio, weight) ;
      // 6j80
      if(njet80>=6)
        _etmiss_HT_6j80->fill( ratio, weight) ;
      // 7j80
      if(njet80>=7)
        _etmiss_HT_7j80->fill( ratio, weight) ;
      // 8j80
      if(njet80>=8)
        _etmiss_HT_8j80->fill( ratio, weight) ;

    }

    //@}

    void finalize() {
      double norm = crossSection()/femtobarn*4.7/sumOfWeights();

      scale(_etmiss_HT_7j55,2.*norm);
      scale(_etmiss_HT_8j55,2.*norm);
      scale(_etmiss_HT_9j55,2.*norm);
      scale(_etmiss_HT_6j80,2.*norm);
      scale(_etmiss_HT_7j80,2.*norm);
      scale(_etmiss_HT_8j80,2.*norm);

      scale(_hist_njet55,norm);
      scale(_hist_njet80,norm);

      scale(_count_7j55,norm);
      scale(_count_8j55,norm);
      scale(_count_9j55,norm);
      scale(_count_6j80,norm);
      scale(_count_7j80,norm);
      scale(_count_8j80,norm);
    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _etmiss_HT_7j55;
    Histo1DPtr _etmiss_HT_8j55;
    Histo1DPtr _etmiss_HT_9j55;
    Histo1DPtr _etmiss_HT_6j80;
    Histo1DPtr _etmiss_HT_7j80;
    Histo1DPtr _etmiss_HT_8j80;

    Histo1DPtr _hist_njet55;
    Histo1DPtr _hist_njet80;

    Histo1DPtr _count_7j55;
    Histo1DPtr _count_8j55;
    Histo1DPtr _count_9j55;
    Histo1DPtr _count_6j80;
    Histo1DPtr _count_7j80;
    Histo1DPtr _count_8j80;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1117704);

}
