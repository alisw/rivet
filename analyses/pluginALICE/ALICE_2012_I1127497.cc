// -*- C++ -*-
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/AliceCommon.hh"
#include "Rivet/Projections/AliceCommon.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {

  /// Analysis
  class ALICE_2012_I1127497 : public Analysis {
  
  public:
    
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I1127497);
    
    /// Book histograms and initialise projections before the run
    void init() {
      
      // Declare centrality projection
      if (getOption("cent") == "REF") {
	_mode = 1;
	declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M");
      }
      else if (getOption("cent") == "IMP") {
	_mode = 2;
	declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality", "V0M", "V0M_IMP");
      }
      
      // Charged final states with |eta| < 0.5 and pT > 150 MeV
      const Cut& cut = Cuts::abseta < 0.5 && Cuts::pT > 150*MeV;
      const ChargedFinalState cfs(cut);
      addProjection(cfs,"CFS");
      
      // Loop over all histograms
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
	
	// Initialize PbPb objects
	_histNch[PBPB][ihist] = bookHisto1D(ihist+1, 1, 1);
	
	std::string nameCounterPbPb = "Counter_PbPb_" + std::to_string(ihist);
	_counterSOW[PBPB][ihist] = bookCounter(nameCounterPbPb, "Sum of weights counter for PbPb");
	
	std::string nameCounterNcoll = "Counter_Ncoll_" + std::to_string(ihist);
	_counterNcoll[ihist] = bookCounter(nameCounterNcoll, "Ncoll counter for PbPb");
	
	// Initialize pp objects. In principle, only one pp histogram would be needed since 
	// centrality does not make any difference here. However, in some cases in this analysis 
	// the binning differ from each other, so this is easy-to-implement way to account for that.
	std::string namePP = _histNch[PBPB][ihist]->name() + "-pp";
	_histNch[PP][ihist] = bookHisto1D(namePP, refData(ihist+1, 1, 1)); // binning taken from ref data
	
	std::string nameCounterpp = "Counter_pp_" + std::to_string(ihist);
	_counterSOW[PP][ihist] = bookCounter(nameCounterpp, "Sum of weights counter for pp");
	
      }
      
      // Centrality regions keeping boundaries for a certain region. 
      // Note, that some regions overlap with other regions.
      _centrRegions.clear();
      _centrRegions = {{0., 5.},   {5., 10.},  {10., 20.},
		       {20., 30.}, {30., 40.}, {40., 50.},
		       {50., 60.}, {60., 70.}, {70., 80.},
		       {0., 10.},  {0., 20.},  {20., 40.},
		       {40., 60.}, {40., 80.}, {60., 80.}};
      
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const double weight = event.weight();
      
      // Final state particles with at least pT = 150 MeV in eta range of |eta| < 0.5
      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      Particles chargedParticles = charged.particlesByPt();
      
      // Check type of event. This may be not the perfect way to check for the type of event as 
      // there might be some weird conditions hidden inside. For example some HepMC versions check 
      // if number of hard collisions is equal to 0 and assign 'false' in that case, which is usually wrong.
      // This might be changed in the future
      if (event.genEvent()->heavy_ion()) {
	
	// Prepare centrality projection and value
	const CentralityProjection& centrProj = apply<CentralityProjection>(event, (_mode == 1 ? "V0M" : "V0M_IMP"));
	double centr = centrProj();
	// Veto event for too large centralities since those are not used in the analysis at all
	if ((centr < 0.) || (centr > 80.))
	  vetoEvent;
	
	// Fill the right PbPb histograms and add weights based on centrality value
	for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
          if (inRange(centr, _centrRegions[ihist].first, _centrRegions[ihist].second)) {
	    _counterSOW[PBPB][ihist]->fill(weight);
	    _counterNcoll[ihist]->fill(event.genEvent()->heavy_ion()->Ncoll(), weight);
	    foreach (const Particle& p, chargedParticles) {
	      float pT = p.pT()/GeV;
	      if (pT < 50.) {
		double pTAtBinCenter = _histNch[PBPB][ihist]->binAt(pT).xMid();
		_histNch[PBPB][ihist]->fill(pT, weight/pTAtBinCenter);
	      }
	    }
	  }
	}
	
      }
      else {
	
	// Fill all pp histograms and add weights
	for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
	  _counterSOW[PP][ihist]->fill(weight);
	  foreach (const Particle& p, chargedParticles) {
	    float pT = p.pT()/GeV;
	    if (pT < 50.) {
	      double pTAtBinCenter = _histNch[PP][ihist]->binAt(pT).xMid();
	      _histNch[PP][ihist]->fill(pT, weight/pTAtBinCenter);
	    }
	  }
	}
	
      }
      
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      
      // Right scaling of the histograms with their individual weights.
      for (size_t itype = 0; itype < EVENT_TYPES; ++itype ) {
	for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
	  if (_counterSOW[itype][ihist]->sumW() > 0.) {
	    scale(_histNch[itype][ihist], (1./_counterSOW[itype][ihist]->sumW() / 2. / M_PI));
	  }
	}
      }
	
      // Postprocessing of the histograms
      for (size_t ihist = 0; ihist < NHISTOS; ++ihist) {
	// If there are entires in histograms for both beam types
	if (_histNch[PP][ihist]->numEntries() > 0 && _histNch[PBPB][ihist]->numEntries() > 0) {
	  // Initialize and fill R_AA histograms
	  _histRAA[ihist] = bookScatter2D(ihist+16, 1, 1);
	  divide(_histNch[PBPB][ihist], _histNch[PP][ihist], _histRAA[ihist]);
	  // Scale by Ncoll. Unfortunately some generators does not provide Ncoll (eg. JEWEL), 
	  // so the following scaling will be done only if there are entries in the counters
	  if (_counterNcoll[ihist]->sumW() > 1e-6 && _counterSOW[PBPB][ihist]->sumW() > 1e-6) {
	    _histRAA[ihist]->scaleY(1. / (_counterNcoll[ihist]->sumW() / _counterSOW[PBPB][ihist]->sumW()));
	  }
	}
      }
      
    }
    
  private:
    
    static const int NHISTOS = 15;
    static const int EVENT_TYPES = 2;
    static const int PP = 0;
    static const int PBPB = 1;
    
    Histo1DPtr _histNch[EVENT_TYPES][NHISTOS];
    CounterPtr _counterSOW[EVENT_TYPES][NHISTOS];
    CounterPtr _counterNcoll[NHISTOS];
    
    Scatter2DPtr _histRAA[NHISTOS];
    std::vector<std::pair<double, double>> _centrRegions;
    
    int _mode;
    
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I1127497);


}
