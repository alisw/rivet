// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {


  /// @brief ATLAS W+b measurement
  class ATLAS_2013_I1219109: public Analysis {
  public:

    ///@brief: Electroweak Wjj production at 8 TeV
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2013_I1219109);
    //@}

    void init() {

      // Get options from the new option system
      _mode = 0;
      if ( getOption("LMODE") == "EL" ) _mode = 1;
      if ( getOption("LMODE") == "MU" ) _mode = 2;

      const FinalState fs;

      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT >= 25*GeV;

      // W finder for electrons and muons
      WFinder wf_mu(fs, cuts, PID::MUON, 0.0*GeV, DBL_MAX, 0.0, 0.1,
                 WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      WFinder wf_el(fs, cuts, PID::ELECTRON, 0.0*GeV, DBL_MAX, 0.0, 0.1,
                 WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
      declare(wf_mu, "WFmu");
      declare(wf_el, "WFel");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(wf_el);
      jet_fs.addVetoOnThisFinalState(wf_mu);
      FastJets fj(jet_fs, FastJets::ANTIKT, 0.4);
      fj.useInvisibles();
      declare(fj, "Jets");
      declare(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "BHadrons");


      // book histograms
      book(_njet     ,1, 1, 1); // dSigma / dNjet
      book(_jet1_bPt ,3, 1, 1); // dSigma / dBjetPt for Njet = 1
      book(_jet2_bPt ,8, 1, 1); // dSigma / dBjetPt for Njet = 2

    }


    void analyze(const Event& event) {

      //  retrieve W boson candidate
      const WFinder& wf_mu = apply<WFinder>(event, "WFmu");
      const WFinder& wf_el = apply<WFinder>(event, "WFel");

      size_t nWmu = wf_mu.size();
      size_t nWel = wf_el.size();

      if (_mode == 0 && !((nWmu == 1 && !nWel) || (!nWmu && nWel == 1)))  vetoEvent; // one W->munu OR W->elnu candidate, otherwise veto
      if (_mode == 1 && !(!nWmu && nWel == 1))  vetoEvent; // one W->elnu candidate, otherwise veto
      if (_mode == 2 && !(nWmu == 1 && !nWel))  vetoEvent; // one W->munu candidate, otherwise veto


      if (   (nWmu? wf_mu : wf_el).bosons().size() != 1 )  vetoEvent; // only one W boson candidate
      if ( !((nWmu? wf_mu : wf_el).mT() > 60.0*GeV) )      vetoEvent;
      //const Particle& Wboson  = wf.boson();


      // retrieve constituent neutrino
      const Particle& neutrino = (nWmu? wf_mu : wf_el).constituentNeutrino();
      if( !(neutrino.pT() > 25*GeV) )  vetoEvent;

      // retrieve constituent lepton
      const Particle& lepton = (nWmu? wf_mu : wf_el).constituentLepton();

      // count good jets, check if good jet contains B hadron
      const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();
      const Jets& jets = apply<JetAlg>(event, "Jets").jetsByPt(25*GeV);
      int goodjets = 0, bjets = 0;
      double bPt = 0.;
      for(const Jet& j : jets) {
        if( (j.abseta() < 2.1) && (deltaR(lepton, j) > 0.5) ) {
          // this jet passes the selection!
          ++goodjets;
          // j.bTagged() uses ghost association which is
          // more elegant, but not what has been used in
          // this analysis originally, will match B had-
          // rons in eta-phi space instead
          for(const Particle& b : bHadrons) {
            if( deltaR(j, b) < 0.3 ) {
              // jet matched to B hadron!
              if(!bPt)  bPt = j.pT() * GeV; // leading b-jet pT
              ++bjets; // count number of b-jets
              break;
            }
          }
        }
      }
      if( goodjets > 2 )  vetoEvent; // at most two jets
      if( !bjets )  vetoEvent; // at least one of them b-tagged

      double njets = double(goodjets);
      double ncomb = 3.0;
      _njet->fill(njets);
      _njet->fill(ncomb);

      if(     goodjets == 1)  _jet1_bPt->fill(bPt);
      else if(goodjets == 2)  _jet2_bPt->fill(bPt);
    }


    void finalize() {
      const double sf = _mode? 1.0 : 0.5;
      const double xs_pb = sf * crossSection() / picobarn  / sumOfWeights();
      const double xs_fb = sf * crossSection() / femtobarn / sumOfWeights();
      scale(_njet,     xs_pb);
      scale(_jet1_bPt, xs_fb);
      scale(_jet2_bPt, xs_fb);
    }

  protected:

    size_t _mode;

  private:

    Histo1DPtr _njet;
    Histo1DPtr _jet1_bPt;
    Histo1DPtr _jet2_bPt;

    //bool _isMuon;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2013_I1219109);

}
