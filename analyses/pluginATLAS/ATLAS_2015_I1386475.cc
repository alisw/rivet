// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/AtlasCommon.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2015_I1386475 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2015_I1386475);
 /// Book histograms and initialise projections before the run
  void init() {

    // The centrality projection.
    declareCentrality(ATLAS::SumET_PB_Centrality(),
                      "ATLAS_pPb_Calib", "SumETPb", "CENT");

    // The trigger projection.
    declare(ATLAS::MinBiasTrigger(), "Trigger");

    // The particles to be analysed.
    declare(ChargedFinalState(Cuts::eta > -2.7 && Cuts::eta < 2.7 &&
                              Cuts::pT > 0.1*GeV), "CFS");
    
    // The centrality bins' upper edges.
    centralityBins = {90., 60., 40., 30., 20., 10., 5., 1.};
    for (int i = 0; i < 8; ++i) {
      book(histEta[centralityBins[i]], 2, 1, i + 1);
      book(sow[centralityBins[i]], "_sow"+to_string(i) + "Counter");
    }
  }

  /// Perform the per-event analysis
  void analyze(const Event& event) {
    
    // Apply event triggers.
    if ( !apply<ATLAS::MinBiasTrigger>(event, "Trigger")() ) vetoEvent;

    // We must have direct acces to the centrality projection.
    const CentralityProjection& cent = 
      apply<CentralityProjection>(event,"CENT");
    double c = cent();
    // Find the correct centrality histogram
    auto hItr = histEta.upper_bound(c);
    if (hItr == histEta.end()) return;
    // Find the correct sow.
    auto sItr = sow.upper_bound(c);
    if (sItr == sow.end()) return;
    sItr->second->fill();
    for ( const auto &p : apply<ChargedFinalState>(event,"CFS").particles() )
      hItr->second->fill(p.eta());
  }
    
  /// Finalize
  void finalize() {

    // Scale by the inverse sum of event weights in each centrality
    // bin.
    for (int i = 0; i < 8; ++i)
      histEta[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());

  }

private:

  /// The histograms binned in centrality.
  vector<double> centralityBins;
  map<double,Histo1DPtr> histEta;
  map<double, CounterPtr> sow;





  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2015_I1386475);


}
