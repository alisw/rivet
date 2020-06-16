// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /// @brief CDF Run II analysis: jet \f$ p_T \f$ and \f$ \eta \f$
  ///   distributions in Z + (b) jet production
  /// @author Lars Sonnenschein
  ///
  /// This CDF analysis provides \f$ p_T \f$ and \f$ \eta \f$ distributions of
  /// jets in Z + (b) jet production, before and after tagging.
  class CDF_2006_S6653332 : public Analysis {
  public:

    /// Constructor
    CDF_2006_S6653332()
      : Analysis("CDF_2006_S6653332"),
        _Rjet(0.7), _JetPtCut(20.), _JetEtaCut(1.5), _Lep1PtCut(18.), _Lep2PtCut(10.), _LepEtaCut(1.1)
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      const FinalState fs((Cuts::etaIn(-3.6, 3.6)));
      declare(fs, "FS");

      // Create a final state with any e+e- or mu+mu- pair with
      // invariant mass 76 -> 106 GeV and ET > 20 (Z decay products)
      vector<pair<PdgId,PdgId> > vids;
      vids.push_back(make_pair(PID::ELECTRON, PID::POSITRON));
      vids.push_back(make_pair(PID::MUON, PID::ANTIMUON));
      FinalState fs2((Cuts::etaIn(-3.6, 3.6)));
      InvMassFinalState invfs(fs2, vids, 66*GeV, 116*GeV);
      declare(invfs, "INVFS");

      // Make a final state without the Z decay products for jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(invfs);
      declare(vfs, "VFS");
      declare(FastJets(vfs, FastJets::CDFMIDPOINT, 0.7), "Jets");

      // Book histograms
      book(_sigmaBJet ,1, 1, 1);
      book(_ratioBJetToZ ,2, 1, 1);
      book(_ratioBJetToJet ,3, 1, 1);

     
      book(_sumWeightsWithZ, "sumWeightsWithZ");
      book(_sumWeightsWithZJet, "sumWeightsWithZJet");
    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Check we have an l+l- pair that passes the kinematic cuts
      // Get the Z decay products (mu+mu- or e+e- pair)
      const InvMassFinalState& invMassFinalState = apply<InvMassFinalState>(event, "INVFS");
      const Particles&  ZDecayProducts =  invMassFinalState.particles();

      // Make sure we have at least 2 Z decay products (mumu or ee)
      if (ZDecayProducts.size() < 2) vetoEvent;
      //
      double Lep1Pt = ZDecayProducts[0].pT();
      double Lep2Pt = ZDecayProducts[1].pT();
      double Lep1Eta = ZDecayProducts[0].absrap(); ///< @todo This is y... should be abseta()?
      double Lep2Eta = ZDecayProducts[1].absrap(); ///< @todo This is y... should be abseta()?

      if (Lep1Eta > _LepEtaCut && Lep2Eta > _LepEtaCut) vetoEvent;
      if (ZDecayProducts[0].abspid()==13 && Lep1Eta > 1. && Lep2Eta > 1.) vetoEvent;
      if (Lep1Pt < _Lep1PtCut && Lep2Pt < _Lep2PtCut) vetoEvent;

      _sumWeightsWithZ->fill();

      /// @todo Write out a warning if there are more than two decay products
      FourMomentum Zmom = ZDecayProducts[0].momentum() +  ZDecayProducts[1].momentum();

      // Put all b-quarks in a vector
      /// @todo Use jet contents rather than accessing quarks directly
      Particles bquarks;
      /// @todo is this HepMC wrangling necessary?
      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())){
        if ( std::abs(p->pdg_id()) == PID::BQUARK ) {
          bquarks.push_back(Particle(*p));
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
          ++numJet;
          // Does the jet contain a b-quark?
          /// @todo Use jet contents rather than accessing quarks directly
          bool bjet = false;
          for (const Particle& bquark : bquarks) {
            if (deltaR(jt->rapidity(), jt->phi(), bquark.rapidity(), bquark.phi()) <= _Rjet) {
              bjet = true;
              break;
            }
          } // end loop around b-jets
          if (bjet) {
            numBJet++;
          }
        }
      } // end loop around jets

      if (numJet > 0) _sumWeightsWithZJet->fill();
      if (numBJet > 0) {
        _sigmaBJet->fill(1960.0);
        _ratioBJetToZ->fill(1960.0);
        _ratioBJetToJet->fill(1960.0);
      }

    }


    /// Finalize
    void finalize() {
      MSG_DEBUG("Total sum of weights = " << sumOfWeights());
      MSG_DEBUG("Sum of weights for Z production in mass range = " << dbl(*_sumWeightsWithZ));
      MSG_DEBUG("Sum of weights for Z+jet production in mass range = " << dbl(*_sumWeightsWithZJet));

      scale(_sigmaBJet, crossSection()/sumOfWeights());
      scale(_ratioBJetToZ, 1.0/ *_sumWeightsWithZ);
      scale(_ratioBJetToJet, 1.0/ *_sumWeightsWithZJet);
    }

    //@}


  private:

    /// @name Cuts and counters
    //@{
    double _Rjet;
    double _JetPtCut;
    double _JetEtaCut;
    double _Lep1PtCut;
    double _Lep2PtCut;
    double _LepEtaCut;

    CounterPtr _sumWeightsWithZ;
    CounterPtr _sumWeightsWithZJet;

    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _sigmaBJet;
    Histo1DPtr _ratioBJetToZ;
    Histo1DPtr _ratioBJetToJet;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2006_S6653332);

}
