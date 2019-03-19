// -*- C++ -*-
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {

  class ALICE_2014_I1243865 : public Analysis {
  // @brief Multi Strange Baryon production at mid rapidity in
  // 2.76 TeV Pb--Pb collisions for different centralities.
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2014_I1243865);
    
    // Book histograms and projections etc.
    void init() {
      // Particles of interest.
      declare(ALICE::PrimaryParticles(Cuts::absrap < 0.5),"CFS"); 

      // The event trigger.
      declare(ALICE::V0AndTrigger(), "V0-AND");

      // The centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(),
           "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Xi Baryons.
      for (string str : {"d01-","d02-","d03-","d04-","d05-"}){
        _histPtXi.push_back(bookHisto1D(str+"x01-y01"));
	_histPtXi_bar.push_back(bookHisto1D(str+"x01-y02"));
      }
      
      // Omega Baryons.
      for (string str : {"d06-","d07-","d08-","d09-","d10-"}){
        _histPtOmega.push_back(bookHisto1D(str+"x01-y01"));
        _histPtOmega_bar.push_back(bookHisto1D(str+"x01-y02"));
      }

      // Sum of weights for the centrality intervals.
      for (int i = 0, N = _histPtOmega.size(); i < N; ++i) {
        sow.push_back(bookCounter("sow_" + toString(i)));
      }

      _histXitoPi = bookProfile1D("d14-x01-y01");
      _histOmegatoPi = bookProfile1D("d14-x01-y02");
    }

    void analyze(const Event& event) {
      // Event weight.
      const double weight = event.weight();
      // Event trigger.
      if (!apply<ALICE::V0AndTrigger>(event, "V0-AND")() ) vetoEvent;
      
      // Centrality. 
      const CentralityProjection& cent = apply<CentralityProjection>(event,"V0M");
      const double c = cent();
		
      int centralityclass = -1;
      if(c > 0. && c <= 10) centralityclass = 0;
      if(c > 10. && c <= 20) centralityclass = 1;
      if(c > 20. && c <= 40) centralityclass = 2;
      if(c > 40. && c <= 60) centralityclass = 3;
      if(c > 60. && c <= 80) centralityclass = 4;
      if (centralityclass == -1) vetoEvent;
      // Fill sum of weights
      sow[centralityclass]->fill(weight);
      int nPions = 0;
      int nXi = 0;
      int nOmega = 0;
      for (const auto& p : 
        apply<ALICE::PrimaryParticles>(event,"CFS").particles()) {
        const double pT = p.pT() / GeV;
        switch (p.pid()){
	  case 211:
	    nPions++;
	    break;
           case 3312:
	     _histPtXi[centralityclass]->fill(pT,weight);
	     nXi++;
	     break;
	   case -3312:
	     _histPtXi_bar[centralityclass]->fill(pT,weight);
	     nXi++;
	     break;
	   case 3334:
	     _histPtOmega[centralityclass]->fill(pT,weight);
	     nOmega++;
	   break;
	     case -3334:
	     _histPtOmega_bar[centralityclass]->fill(pT,weight);
	     nOmega++;
	   break;
         }
       }
       // Extract Npart form GenEvent. TODO: Unclear how to do
       // this in HepMC3
       const HepMC::HeavyIon* hi = event.genEvent()->heavy_ion();
       if (hi && nPions != 0){
	 const double npart = hi->Npart_proj() + hi->Npart_targ();
         if (nXi != 0) 
           _histXitoPi->fill(npart, double(nXi) / double(nPions), weight);
         if (nOmega != 0) 
	   _histOmegatoPi->fill(npart, double(nOmega) / double(nPions),
	      weight);
	}
     }

    void finalize() {
      
      for (int i = 0, N = _histPtOmega.size(); i < N; ++i) {
        const double s = 1./sow[i]->sumW();
        _histPtXi[i]->scaleW(s);
        _histPtXi_bar[i]->scaleW(s);
        _histPtOmega[i]->scaleW(s);
        _histPtOmega_bar[i]->scaleW(s);
      }
   }

  private:
    vector<Histo1DPtr> _histPtXi;
    vector<Histo1DPtr> _histPtXi_bar;
    vector<Histo1DPtr> _histPtOmega;
    vector<Histo1DPtr> _histPtOmega_bar;
    vector<CounterPtr> sow;
    Profile1DPtr _histXitoPi;
    Profile1DPtr _histOmegatoPi;
};

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2014_I1243865);
}

