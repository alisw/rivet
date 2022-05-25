// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
namespace Rivet {


  /// @brief ALICE PbPb at 5.02 TeV eta distributions.
  class ALICE_2016_I1507090 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2016_I1507090);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      // Centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality",
        "V0M","V0M");
      // Projections for the 2-out-of-3 trigger.
      declare(ChargedFinalState( (Cuts::eta > 2.8 && Cuts::eta < 5.1) &&
         Cuts::pT > 0.1*GeV), "VZERO1");
      declare(ChargedFinalState( (Cuts::eta > -3.7 && Cuts::eta < -1.7) &&
	 Cuts::pT > 0.1*GeV), "VZERO2");
      declare(ChargedFinalState(Cuts::abseta < 1. && Cuts::pT > 0.15*GeV),
        "SPD");

      // Primary particles.
      declare(ALICE::PrimaryParticles(Cuts::abseta < 5.6),"APRIM");

      // The centrality bins and the corresponding histograms.
      centralityBins =
      { 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 };
      for (int i = 0, N = centralityBins.size(); i < N; ++i) {
        book(histEta[centralityBins[i]], 1, 1, i + 1);
        book(sow[centralityBins[i]], "sow_" + toString(i));
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Trigger projections.
      const ChargedFinalState& vz1 =
        applyProjection<ChargedFinalState>(event,"VZERO1");
      const ChargedFinalState& vz2 =
        applyProjection<ChargedFinalState>(event,"VZERO2");
      const ChargedFinalState& spd =
        applyProjection<ChargedFinalState>(event,"SPD");
      int fwdTrig = (vz1.particles().size() > 0 ? 1 : 0);
      int bwdTrig = (vz2.particles().size() > 0 ? 1 : 0);
      int cTrig = (spd.particles().size() > 0 ? 1 : 0);

      if (fwdTrig + bwdTrig + cTrig < 2) vetoEvent;
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      double c = cent();
      // Find the correct centrality histogram
      auto hItr = histEta.upper_bound(c);
      if (hItr == histEta.end()) return;
      // Find the correct sow.
      auto sItr = sow.upper_bound(c);
      if (sItr == sow.end()) return;
      sItr->second->fill();

      // Fill the histograms.
      for ( const auto& p :
        applyProjection<ALICE::PrimaryParticles>(event,"APRIM").particles() )
	if(p.abscharge() > 0) hItr->second->fill(p.eta());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int i = 0, N = centralityBins.size(); i < N; ++i)
        histEta[centralityBins[i]]->scaleW(1./sow[centralityBins[i]]->sumW());

    }

    //@}


    /// @name Histograms
    //@{
    vector<double> centralityBins;
    map<double,Histo1DPtr> histEta;
    map<double, CounterPtr> sow;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALICE_2016_I1507090);


}
