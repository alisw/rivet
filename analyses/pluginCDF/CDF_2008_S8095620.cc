// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"

namespace Rivet {


  /// @brief CDF Run II Z + b-jet cross-section measurement
  class CDF_2008_S8095620 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2008_S8095620);


    /// @name Analysis methods
    //@{

    void init() {
      // Set up projections
      const FinalState fs((Cuts::etaIn(-3.2, 3.2)));
      declare(fs, "FS");
      // Create a final state with any e+e- or mu+mu- pair with
      // invariant mass 76 -> 106 GeV and ET > 18 (Z decay products)
      vector<pair<PdgId,PdgId> > vids;
      vids.push_back(make_pair(PID::ELECTRON, PID::POSITRON));
      vids.push_back(make_pair(PID::MUON, PID::ANTIMUON));
      FinalState fs2((Cuts::etaIn(-3.2, 3.2)));
      InvMassFinalState invfs(fs2, vids, 76*GeV, 106*GeV);
      declare(invfs, "INVFS");
      // Make a final state without the Z decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      declare(vfs, "VFS");
      declare(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");

      // Book histograms
      book(_dStot    ,1, 1, 1);
      book(_dSdET    ,2, 1, 1);
      book(_dSdETA   ,3, 1, 1);
      book(_dSdZpT   ,4, 1, 1);
      book(_dSdNJet  ,5, 1, 1);
      book(_dSdNbJet ,6, 1, 1);

      book(_sumWeightSelected,"sumWeightSelected");
    }


    // Do the analysis
    void analyze(const Event& event) {
      // Check we have an l+l- pair that passes the kinematic cuts
      // Get the Z decay products (mu+mu- or e+e- pair)
      const InvMassFinalState& invMassFinalState = apply<InvMassFinalState>(event, "INVFS");
      const Particles&  ZDecayProducts =  invMassFinalState.particles();

      // make sure we have 2 Z decay products (mumu or ee)
      if (ZDecayProducts.size() < 2) vetoEvent;
      //new cuts
      double Lep1Pt = ZDecayProducts[0].perp();
      double Lep2Pt = ZDecayProducts[1].perp();
      double Lep1Eta = fabs(ZDecayProducts[0].rapidity());
      double Lep2Eta = fabs(ZDecayProducts[1].rapidity());

      if (Lep1Eta > _LepEtaCut || Lep2Eta > _LepEtaCut) vetoEvent;

      if (ZDecayProducts[0].abspid()==13 &&
          ((Lep1Eta > 1.5 || Lep2Eta > 1.5) || (Lep1Eta > 1.0 && Lep2Eta > 1.0))) {
        vetoEvent;
      }

      if (Lep1Pt > Lep2Pt) {
        if (Lep1Pt < _Lep1PtCut || Lep2Pt < _Lep2PtCut) vetoEvent;
      }
      else {
        if (Lep1Pt < _Lep2PtCut || Lep2Pt < _Lep1PtCut) vetoEvent;
      }

      _sumWeightSelected->fill();
      /// @todo: write out a warning if there are more than two decay products
      FourMomentum Zmom = ZDecayProducts[0].momentum() +  ZDecayProducts[1].momentum();

      // Put all b-quarks in a vector
      /// @todo Use a b-hadron search rather than b-quarks for tagging
      Particles bquarks;
      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {
        if (std::abs(p->pdg_id()) == PID::BQUARK) {
          bquarks += Particle(*p);
        }
      }

      // Get jets
      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      MSG_DEBUG("Jet multiplicity before any pT cut = " << jetpro.size());

      const PseudoJets& jets = jetpro.pseudoJetsByPt();
      MSG_DEBUG("jetlist size = " << jets.size());

      int numBJet = 0;
      int numJet  = 0;
      // for each b-jet plot the ET and the eta of the jet, normalise to the total cross section at the end
      // for each event plot N jet and pT(Z), normalise to the total cross section at the end
      for (PseudoJets::const_iterator jt = jets.begin(); jt != jets.end(); ++jt) {
        // select jets that pass the kinematic cuts
        if (jt->perp() > _JetPtCut && fabs(jt->rapidity()) <= _JetEtaCut) {
          numJet++;
          // does the jet contain a b-quark?
          bool bjet = false;
          for (const Particle& bquark :  bquarks) {
            if (deltaR(jt->rapidity(), jt->phi(), bquark.rapidity(),bquark.phi()) <= _Rjet) {
              bjet = true;
              break;
            }
          } // end loop around b-jets
          if (bjet) {
            numBJet++;
            _dSdET->fill(jt->perp());
            _dSdETA->fill(fabs(jt->rapidity()));
          }
        }
      } // end loop around jets

      // wasn't asking for b-jets before!!!!
      if(numJet > 0 && numBJet > 0) _dSdNJet->fill(numJet);
      if(numBJet > 0) {
        _dStot->fill(1960.0);
        _dSdNbJet->fill(numBJet);
        _dSdZpT->fill(Zmom.pT());
      }
    }


    // Finalize
    void finalize() {
      // normalise histograms
      // scale by 1 / the sum-of-weights of events that pass the Z cuts
      // since the cross sections are normalized to the inclusive
      // Z cross sections.
      double Scale = 1.0;
      if (_sumWeightSelected->val() != 0.0) Scale = 1.0/dbl(*_sumWeightSelected);
      scale(_dStot,Scale);
      scale(_dSdET,Scale);
      scale(_dSdETA,Scale);
      scale(_dSdNJet,Scale);
      scale(_dSdNbJet,Scale);
      scale(_dSdZpT,Scale);
    }

    //@}


  private:

    /// @name Cuts
    /// @{
    const double _Rjet = 0.7;
    const double _JetPtCut = 20;
    const double _JetEtaCut = 1.5;
    const double _Lep1PtCut = 18;
    const double _Lep2PtCut = 10;
    const double _LepEtaCut = 3.2;
    /// @}

    /// Counter
    CounterPtr _sumWeightSelected;

    /// @name Histograms
    /// @{
    Histo1DPtr _dStot;
    Histo1DPtr _dSdET;
    Histo1DPtr _dSdETA;
    Histo1DPtr _dSdNJet;
    Histo1DPtr _dSdNbJet;
    Histo1DPtr _dSdZpT;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2008_S8095620, CDF_2008_I806082);

}
