// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/CentralityProjection.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include "Rivet/Projections/HepMCHeavyIon.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {

  class ALICE_2014_I1243865 : public Analysis {

  // @brief Multi-strange baryon production at mid-rapidity in 2.76 TeV Pb--Pb collisions
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2014_I1243865);

    
    // Book histograms and projections etc.
    void init() {
      // Particles of interest.
      declare(ALICE::PrimaryParticles(Cuts::absrap < 0.5),"CFS");

      // The event trigger.
      declare(ALICE::V0AndTrigger(), "V0-AND");

      // The centrality projection.
      declareCentrality(ALICE::V0MMultiplicity(),
           "ALICE_2015_PBPBCentrality", "V0M", "V0M");

      // Access the HepMC heavy ion info
      declare(HepMCHeavyIon(), "HepMC");

      // Xi Baryons.
      size_t ixi = 0;
      _histPtXi.resize(5);
      _histPtXi_bar.resize(5);
      for (string str : {"d01-","d02-","d03-","d04-","d05-"}){
        book(_histPtXi[ixi], str+"x01-y01");
        book(_histPtXi_bar[ixi], str+"x01-y02");
        ixi += 1;
      }

      // Omega Baryons.
      size_t iom = 0;
      _histPtOmega.resize(5);
      _histPtOmega_bar.resize(5);
      for (string str : {"d06-","d07-","d08-","d09-","d10-"}){
        book(_histPtOmega[iom], str+"x01-y01");
        book(_histPtOmega_bar[iom], str+"x01-y02");
        iom += 1;
      }

      // Sum of weights for the centrality intervals.
      sow.resize(_histPtOmega.size());
      for (int i = 0, N = _histPtOmega.size(); i < N; ++i) {
        book(sow[i], "sow_" + toString(i));
      }

      book(_histXitoPi, "d14-x01-y01");
      book(_histOmegatoPi, "d14-x01-y02");
    }


    void analyze(const Event& event) {
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
      sow[centralityclass]->fill();
      int nPions = 0;
      int nXi = 0;
      int nOmega = 0;
      for (const Particle& p : apply<ALICE::PrimaryParticles>(event,"CFS").particles()) {
        const double pT = p.pT() / GeV;
        switch (p.pid()){
	  case 211:
	    nPions++;
	    break;
           case 3312:
	     _histPtXi[centralityclass]->fill(pT);
	     nXi++;
	     break;
	   case -3312:
	     _histPtXi_bar[centralityclass]->fill(pT);
	     nXi++;
	     break;
	   case 3334:
	     _histPtOmega[centralityclass]->fill(pT);
	     nOmega++;
	   break;
	     case -3334:
	     _histPtOmega_bar[centralityclass]->fill(pT);
	     nOmega++;
	   break;
         }
       }
       // Extract Npart form GenEvent. TODO: Unclear how to do
       // this in HepMC3
      const HepMCHeavyIon & hi = apply<HepMCHeavyIon>(event, "HepMC");
       if ( nPions != 0){
	 const double npart = hi.Npart_proj() + hi.Npart_targ();
         if (nXi != 0) 
           _histXitoPi->fill(npart, double(nXi) / double(nPions));
         if (nOmega != 0) 
	   _histOmegatoPi->fill(npart, double(nOmega) / double(nPions));
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
  RIVET_DECLARE_PLUGIN(ALICE_2014_I1243865);

}

