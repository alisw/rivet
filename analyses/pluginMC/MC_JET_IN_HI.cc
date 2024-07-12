// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "Rivet/Tools/AtlasCommon.hh"

namespace Rivet {

  class MC_JET_IN_HI : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_JET_IN_HI);

    //@}


  public:

    string ts(int in) {
      std::stringstream ss;
      ss << in;
      return ss.str();
    }
    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Declare centrality projection - we use the ATLAS PbPb definition
      // to be able to compare to data.
      declareCentrality(ATLAS::SumET_PBPB_Centrality(),"ATLAS_PBPB_CENTRALITY",
		      "sumETFwd","sumETFwd");
      // The final state where jets are found
      FinalState fs(Cuts::abseta < 2.5);
      declare(fs, "FS");
 
      ZFinder zfinder(fs, Cuts::abseta < 2.5 && Cuts::pT > 30*GeV, PID::MUON, 80*GeV, 100*GeV, 0.2,
                      ZFinder::ChargedLeptons::PROMPT, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
      declare(zfinder, "ZFinder");

      // Z+jet jet collections
      declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.3), "JetsAK3");
      declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.5), "JetsAK5");
      declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.7), "JetsAK7");
      declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.9), "JetsAK9");

      jetFinders = {"JetsAK3", "JetsAK5", "JetsAK7", "JetsAK9"};
      
      h_zpT.resize(jetFinders.size());
      h_jetpT.resize(jetFinders.size());
      for (size_t i = 0; i < jetFinders.size(); ++i) {
        string s = jetFinders[i];
        book(h_zpT[i], s + "zpT",logspace(50, 1.0,1000));
        book(h_jetpT[i], s + "jetpT",logspace(50, 1.0,1000));
      }
      book(incSow, "incSow");

      centData = {0., 0.2, 0.4, 0.6, 0.8,};
      for (size_t i = 0; i < centData.size(); ++i) {
        book(c_jetpT[centData[i]], "cjetpT" + ts(i),logspace(100, 10.0,1000));
        book(c_zpT[centData[i]], "czpt" + ts(i),logspace(100, 10.0,1000));
	      book(sow[centData[i]], "sow_" + ts(i));
      }
    }

    bool isBackToBack_zj(const ZFinder& zf, const fastjet::PseudoJet& psjet) {
      const FourMomentum& z = zf.bosons()[0].momentum();
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(z, jmom) > 7.*M_PI/8. );
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      // Get the Z
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      Particle z = zfinder.bosons()[0];
      Particle l1 = zfinder.constituents()[0];
      Particle l2 = zfinder.constituents()[1];

      // Require a high-pT Z (and constituents)
      if (l1.pT() < 10*GeV || l2.pT() < 10*GeV || z.pT() < 60*GeV) vetoEvent;

      // Get the centrality
      const double c = apply<CentralityProjection>(event,"sumETFwd")();
      auto jetpTItr = c_jetpT.upper_bound(c);
      auto zpTItr = c_zpT.upper_bound(c);
      auto sowItr = sow.upper_bound(c);
      if (jetpTItr == c_jetpT.end() || zpTItr == c_zpT.end() || 
        sowItr == sow.end()) vetoEvent;
      sowItr->second->fill();
      incSow->fill();
      // Get the  jets
      for (size_t i = 0; i < jetFinders.size(); ++i ) {
        const PseudoJets& psjets = apply<FastJets>(event, 
	  jetFinders[i]).pseudoJetsByPt(30.0*GeV);
        if (!psjets.empty()) {
        // Get the leading jet and make sure it's back-to-back with the Z
        const fastjet::PseudoJet& j0 = psjets[0];
	if (isBackToBack_zj(zfinder, j0)) {
	    // Fill the centrality inclusive histograms
	    h_zpT[i]->fill(z.pT());
	    h_jetpT[i]->fill(j0.perp());
	    // Fill centrality dept histograms only for R = 0.3
	    if (i == 0) {
              jetpTItr->second->fill(j0.perp());
	      zpTItr->second->fill(z.pT());
	    }
	  }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
       for(size_t i = 0; i < jetFinders.size(); ++i) {
         h_zpT[i]->scaleW(1./incSow->sumW());
         h_jetpT[i]->scaleW(1./incSow->sumW());
       }
       for (size_t i = 0; i < centData.size(); ++i) {
         c_jetpT[centData[i]]->scaleW(1./sow[centData[i]]->sumW());
         c_zpT[centData[i]]->scaleW(1./sow[centData[i]]->sumW());
       }
    }


    //@}


  private:
    vector<string> jetFinders;
    // Centrality inclusive histograms
    vector<Histo1DPtr> h_zpT;
    vector<Histo1DPtr> h_jetpT;
    CounterPtr incSow;
    // Centrality intervals
    vector<double> centData;
    // Centrality binned histograms
    map<double, Histo1DPtr> c_jetpT;
    map<double, Histo1DPtr> c_zpT;
    map<double, CounterPtr> sow;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_JET_IN_HI);

}
