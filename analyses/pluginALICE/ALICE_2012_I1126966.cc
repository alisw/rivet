//-*- C++ -*-
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {
  
  
  /// Pion, Kaon, and Proton Production in 0-5%
  ///  central Pb--Pb Collisions at 2.76 TeV
  class ALICE_2012_I1126966 : public Analysis {
  public:
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I1126966);
    
    /// Book histograms and initialise projections before the run
    void init() {
      // Particles of interest.
      declare(ALICE::PrimaryParticles(Cuts::absrap < 0.5),"CFS"); 

      // The event trigger.
      declare(ALICE::V0AndTrigger(), "V0-AND");
      // The centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(),
           "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Invariant pT distributions.
      _histPtPi = bookHisto1D("d01-x01-y01"); //pi+
      _histPtPibar = bookHisto1D("d01-x01-y02");// pi-
      _histPtKaon = bookHisto1D("d02-x01-y01"); //K+
      _histPtKaonbar = bookHisto1D("d02-x01-y02"); //K-
      _histPtProton = bookHisto1D("d03-x01-y01"); //P+
      _histPtProtonbar = bookHisto1D("d03-x01-y02"); //P-
      // Yield histograms.
      _histNpi = bookHisto1D("d04-x01-y01");
      _histNpibar = bookHisto1D("d04-x01-y02");
      _histNKaon = bookHisto1D("d04-x01-y03");
      _histNKaonbar = bookHisto1D("d04-x01-y04");
      _histNproton = bookHisto1D("d04-x01-y05");
      _histNprotonbar =bookHisto1D("d04-x01-y06");
      // Sum of weights of triggered events.
      sow = bookCounter("sow");

  }
    
      /// Perform the per-event analysis
    
    void analyze(const Event& event) {
      // Event weight.
      const double weight = event.weight();
      // Analysis only considers 0-5% central events
      if (apply<CentralityProjection>(event,"V0M")() > 5.0)
        vetoEvent;
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;
      
      sow->fill(weight);
      // ID particles counters for this event.
      int Npi=0;
      int Npibar=0;
      int NKaon=0;
      int NKaonbar=0;
      int Nproton=0;
      int Nprotonbar=0;

      for (const auto& p : 
	apply<ALICE::PrimaryParticles>(event,"CFS").particles()) {
          const double pWeight = weight / p.pT() / 2. / M_PI;	      
          switch (p.pid()) {
            case 211: // pi+
	      Npi++;
              _histPtPi->fill(p.pT()/GeV,  pWeight);
              break;
	    case -211: //pi-
	      Npibar++;
	      _histPtPibar->fill(p.pT()/GeV, pWeight);
              break;
            case 2212: // proton
	      Nproton++;
              _histPtProton->fill(p.pT()/GeV, pWeight);
	      break;
	    case -2212: // p-bar
	      Nprotonbar++;
              _histPtProtonbar->fill(p.pT()/GeV, pWeight);
              break;
            case 321: // K+
	      NKaon++;
              _histPtKaon->fill(p.pT()/GeV,  pWeight);
	      break;
	    case -321: // K-
	      NKaonbar++;
              _histPtKaonbar->fill(p.pT()/GeV,  pWeight);
             break;
        }
      } // Particle loop ends.

      // Fill yield histograms.

      _histNpi->fill(0.0, Npi, weight);
      _histNpibar->fill(0.0, Npibar, weight);
      _histNKaon->fill(0.0, NKaon, weight);
      _histNKaonbar->fill(0.0, NKaonbar, weight);
      _histNproton->fill(0.0, Nproton, weight);
      _histNprotonbar->fill(0.0, Nprotonbar, weight);
    }


    void finalize() {
       const double s = 1./sow->sumW();
       _histPtPi->scaleW(s);
       _histPtPibar->scaleW(s);
       _histPtKaon->scaleW(s);
       _histPtKaonbar->scaleW(s);
       _histPtProton->scaleW(s);
       _histPtProtonbar->scaleW(s);
       _histNpi->scaleW(s);
       _histNpibar->scaleW(s);
       _histNKaon->scaleW(s);
       _histNKaonbar->scaleW(s);
       _histNproton->scaleW(s);
       _histNprotonbar->scaleW(s);

}
    
  private:
    
      // pT histograms
    Histo1DPtr _histPtPi;
    Histo1DPtr _histPtKaon;
    Histo1DPtr _histPtProton;
    Histo1DPtr _histPtPibar;
    Histo1DPtr _histPtKaonbar;
    Histo1DPtr _histPtProtonbar;
    Histo1DPtr _histNpi;
    Histo1DPtr _histNpibar;
    Histo1DPtr _histNKaon;
    Histo1DPtr _histNKaonbar;
    Histo1DPtr _histNproton;
    Histo1DPtr _histNprotonbar;
    CounterPtr sow;

  };
  
  
    // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I1126966);
  
  
}
