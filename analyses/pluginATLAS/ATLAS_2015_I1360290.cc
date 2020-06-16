// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AtlasCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2015_I1360290 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1360290);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      // Centrality projection.
      declareCentrality(ATLAS::SumET_PBPB_Centrality(), "ATLAS_PBPB_CENTRALITY",
	"sumETFwd","sumETFwd");
      // Trigger projection.
      declare(ATLAS::MinBiasTrigger(),"Trigger");
      // The measured final state.
      declare(ChargedFinalState (Cuts::abseta < 2. &&
        Cuts::pT > 0.5*GeV && Cuts::pT < 150.0*GeV), "CFS");
      
      taa = {26.3, 20.6, 14.4, 8.73, 5.05, 2.70, 1.34, 0.41};
      centData = {5., 10., 20., 30., 40., 50., 60., 80.};
      
      for (int i = 0, N = centData.size(); i < N; ++i) {
	// eta hists starts from table 55 ( first 1.7 < pT < 2.0)
        book(histEta1[centData[i]], 55 + i, 1, 1);
	// From table 64, 6.7 < pT < 7.7
	book(histEta2[centData[i]], 64 + i, 1, 1 );
	// From table 73, 19.9 < pT < 22.8
	book(histEta3[centData[i]], 73 + i, 1, 1 );
	// From table 82, 59.8 < pT < 94.8
	book(histEta4[centData[i]], 82 + i, 1, 1 );
	// pt hists starts from table 2 on hepmc, |eta| < 2.0
	book(histpT[centData[i]], 2 + i, 1, 1);
	// keep track of sow in centrality bins.
	book(sow[centData[i]], "sow_" + toString(i));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if ( !apply<ATLAS::MinBiasTrigger>(event, "Trigger")() ) vetoEvent;
     const CentralityProjection& cent = apply<CentralityProjection>(event,"sumETFwd");
      double c = cent();
      // Find the correct centrality histograms
      auto hItr1 = histEta1.upper_bound(c);
      if (hItr1 == histEta1.end()) return;
      auto hItr2 = histEta2.upper_bound(c);
      if (hItr2 == histEta2.end()) return;
      auto hItr3 = histEta3.upper_bound(c);
      if (hItr3 == histEta3.end()) return;
      auto hItr4 = histEta4.upper_bound(c);
      if (hItr4 == histEta4.end()) return;
      auto hpTItr = histpT.upper_bound(c);
      if (hpTItr == histpT.end()) return;
      // Find the correct sow.
      auto sItr = sow.upper_bound(c);
      if (sItr == sow.end()) return;
      sItr->second->fill();

      for (const auto& p : apply<ChargedFinalState>(event,"CFS").particles()) {
        double pT = p.pT();
	double eta = abs(p.eta());
	if (pT > 1.7 && pT < 2.0) hItr1->second->fill(eta, 0.5);
	else if (pT > 6.7 && pT < 7.7) hItr2->second->fill(eta, 0.5);
	else if (pT > 19.9 && pT < 22.8) hItr3->second->fill(eta, 0.5);
	else if (pT > 59.8 && pT < 94.8) hItr4->second->fill(eta, 0.5);
	if (eta < 2) hpTItr->second->fill(pT, 1.0/2./M_PI/pT/4.);
      
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int i = 0, N = centData.size(); i < N; ++i) {
        histEta1[centData[i]]->scaleW(1./sow[centData[i]]->sumW());
        histEta2[centData[i]]->scaleW(1./sow[centData[i]]->sumW());
        histEta3[centData[i]]->scaleW(1./sow[centData[i]]->sumW());
        histEta4[centData[i]]->scaleW(1./sow[centData[i]]->sumW());
        histpT[centData[i]]->scaleW(1./sow[centData[i]]->sumW()/taa[i]);

      }

    }

    //@}


    /// @name Histograms
    //@{
    // The centrality binned histograms
    map<double, Histo1DPtr> histEta1;
    map<double, Histo1DPtr> histEta2;
    map<double, Histo1DPtr> histEta3;
    map<double, Histo1DPtr> histEta4;
    map<double, Histo1DPtr> histpT;
    map<double, CounterPtr> sow;


    vector<double> centData;
    vector<double> taa;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1360290);


}
