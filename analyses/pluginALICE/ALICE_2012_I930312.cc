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
  class ALICE_2012_I930312 : public Analysis {
  
  public:
    
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I930312);
    
    /// Book histograms and initialise projections before the run
    void init() {

      // Declare centrality projection
      declareCentrality(ALICE::V0MMultiplicity(), "ALICE_2015_PBPBCentrality",
        "V0M", "V0M");

      // Charged final states with |eta| < 1.0 and 8 < pT < 15 GeV/c for 
      // trigger particles
      const Cut& cutTrigger = 
        Cuts::abseta < 1.0 && Cuts::pT > 8*GeV && Cuts::pT < 15*GeV;
      const ChargedFinalState cfsTrigger(cutTrigger);
      addProjection(cfsTrigger,"CFSTrigger");
      
      // Set limit values of pT bins
      pt_limits[0] = 3.;
      pt_limits[1] = 4.;
      pt_limits[2] = 6.;
      pt_limits[3] = 8.;
      pt_limits[4] = 10.;
      
      // Charged final states with |eta| < 1.0 and different pT bins for 
      // associated particles.
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	Cut mycut = Cuts::abseta < 1.0 && Cuts::pT > pt_limits[ipt]*GeV && 
	  Cuts::pT < pt_limits[ipt + 1]*GeV;
	declare(ChargedFinalState(mycut), "CFSAssoc" + std::to_string(ipt));
      }
      
      // Create event strings
      event_string[0] = "pp";
      event_string[1] = "central";
      event_string[2] = "peripheral";
      event_string[3] = "other";
      
      // For each event type
      for (int itype = 0; itype < EVENT_TYPES; itype++) {
	// For each pT range
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  // Initialize yield histograms
	  _histYield[itype][ipt] = bookHisto1D("Yield_" + event_string[itype]
	    + "_" + std::to_string(ipt), 36, -0.5 * M_PI, 1.5 * M_PI, 
	    "Associated particle per trigger particle yield", 
	    "$\\Delta\\eta$ (rad)", "$1 / N_{trig} dN_{assoc} / d\\Delta\\eta$ (rad$^-1$)");
	 _histYieldBkgRemoved[itype][ipt] = bookHisto1D("Yield_" + 
	   event_string[itype] + "_nobkg_" + std::to_string(ipt), 36, -0.5*M_PI,
	   1.5 * M_PI, "Associated particle per trigger particle yield no bkg", 
	   "$\\Delta\\eta$ (rad)", "$1 / N_{trig} dN_{assoc} / d\\Delta\\eta$ (rad$^-1$)");
	 }
      }
      
      // Histogram for counting trigger particles for each event type
      _histTriggerCounter = bookHisto1D("Trigger", EVENT_TYPES, 0.0, 
        EVENT_TYPES, "Trigger counter", "event type", "N");
      
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      
      const double weight = event.weight();
      
      // Create charged final state for trigger particle
      const ChargedFinalState& triggerFinalState = 
        applyProjection<ChargedFinalState>(event, "CFSTrigger");
      Particles triggerParticles = triggerFinalState.particlesByPt();

      // Create charged final state for associated particle      
      ChargedFinalState associatedFinalState[PT_BINS];
      Particles associatedParticles[PT_BINS];
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	associatedFinalState[ipt] = applyProjection<ChargedFinalState>(event, 
	  "CFSAssoc" + std::to_string(ipt));
	associatedParticles[ipt] = associatedFinalState[ipt].particlesByPt();
      }
      
      // Check event type
      if (event.genEvent()->heavy_ion()) {
	// Prepare centrality projection and value
	const CentralityProjection& centrProj = 
          apply<CentralityProjection>(event, "V0M");
	double centr = centrProj();
	
	// Set the flag for the type of the event
	if (centr > 0.0 && centr < 5.0)
	  event_type = 1; // PbPb, central
	else if (centr > 60.0 && centr < 90.0)
	  event_type = 2; // PbPb, peripherial
	else
	  event_type = 3; // PbPb, other
      }
      else {
	event_type = 0; // pp
      }
      
      // Veto event if not valid event type
      if (event_type == 3)
	vetoEvent;
      
      // Fill trigger histogram for a proper event type
      _histTriggerCounter->fill(event_type, triggerParticles.size());
      
      // Loop over trigger particles
      for (const auto& triggerParticle : triggerParticles) {
	// For each pt bin
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  // Loop over associated particles
	  for (const auto& associatedParticle : associatedParticles[ipt]) {
	    // If associated and trigger particle are not the same particles.
	    if (associatedParticle != triggerParticle) {
	      // Test trigger particle.
	      if (triggerParticle.pt() > associatedParticle.pt()) {
		// Calculate delta phi in range (-0.5*PI, 1.5*PI).
		double dPhi = triggerParticle.phi() - 
		  associatedParticle.phi();
		while (dPhi > 1.5 * M_PI)  { dPhi -= 2 * M_PI; }
		while (dPhi < -0.5 * M_PI) { dPhi += 2 * M_PI; }
		// Fill yield histogram for calculated delta phi
		_histYield[event_type][ipt]->fill(dPhi, weight);
	      }
	    }
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      
      // Check for the reentrant finalize
      bool pp_available = false, PbPb_available = false;
      reentrant_flag = false;
      // For each event type
      for (int itype = 0; itype < EVENT_TYPES; itype++) {
	// For each pT range
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  if (_histYield[itype][ipt]->numEntries() > 0)
	    itype == 0 ? pp_available = true : PbPb_available = true;
	}
      }
      reentrant_flag = (pp_available && PbPb_available);
      
      // Postprocessing of the histograms
      if (reentrant_flag == true) {
	
        // Initialize IAA and ICP histograms
        _histIAA[0] = bookScatter2D(1, 1, 1);
        _histIAA[1] = bookScatter2D(2, 1, 1);
        _histIAA[2] = bookScatter2D(5, 1, 1);

        _histIAA[3] = bookScatter2D(3, 1, 1);
        _histIAA[4] = bookScatter2D(4, 1, 1);
        _histIAA[5] = bookScatter2D(6, 1, 1);
	// Variables for near and away side peak calculation
	double nearSide[EVENT_TYPES][PT_BINS] = { {0.0} };
	double awaySide[EVENT_TYPES][PT_BINS] = { {0.0} };
      
	// Variables for background error calculation
	double background[EVENT_TYPES][PT_BINS] = { {0.0} };
	double backgroundError[EVENT_TYPES][PT_BINS] = { {0.0} };
      
	// Variables for integration error calculation
	double scalingFactor[EVENT_TYPES] = {0.0};
	double numberOfEntries[EVENT_TYPES][PT_BINS][2] = { { {0.0} } };
	int numberOfBins[EVENT_TYPES][PT_BINS][2] = { { {0} } };
      
	// For each event type
	for (int itype = 0; itype < EVENT_TYPES; itype++) {
	  // For each pT range
	  for (int ipt = 0; ipt < PT_BINS; ipt++) {
	    // Check if histograms are fine
	    if (_histTriggerCounter->numEntries() == 0 || 
	      _histYield[itype][ipt]->numEntries() == 0) {
	      cout << "There are no entries in one of the histograms" << endl;
	      continue;
	    }
	  
	    // Scale yield histogram
	    if ((_histTriggerCounter->bin(itype).sumW() != 0)) {
	      scalingFactor[itype] = 1. / 
	        _histTriggerCounter->bin(itype).sumW();
	      scale(_histYield[itype][ipt], 
	        (1. / _histTriggerCounter->bin(itype).sumW()));
	    }
	  
	    // Calculate background
	    double sum = 0.0;
	    int nbins = 0;
	    for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	      if ((_histYield[itype][ipt]->bin(ibin).xMid() > (-0.5 * M_PI) && 
		   _histYield[itype][ipt]->bin(ibin).xMid() < (-0.5 * M_PI + 0.4)) ||
		  (_histYield[itype][ipt]->bin(ibin).xMid() > (0.5 * M_PI - 0.4) &&
		   _histYield[itype][ipt]->bin(ibin).xMid() < (0.5 * M_PI + 0.4)) ||
		  (_histYield[itype][ipt]->bin(ibin).xMid() > (1.5 * M_PI - 0.4) &&
		   _histYield[itype][ipt]->bin(ibin).xMid() < (1.5 * M_PI))) {
		sum += _histYield[itype][ipt]->bin(ibin).sumW();
		nbins++;
	      }
	    }
	    if (nbins == 0) {
	      std::cout << "Failed to estimate background!" << std::endl;
	      continue;
	    }
	    background[itype][ipt] = sum / nbins;
	  
	    // Calculate background error
	    sum = 0.0;
	    nbins = 0;
	    for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	      if (_histYield[itype][ipt]->bin(ibin).xMid() > (0.5 * M_PI - 0.4) &&
		  _histYield[itype][ipt]->bin(ibin).xMid() < (0.5 * M_PI + 0.4)) {
		sum += (_histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt]) *
		  (_histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt]);
		nbins++;
	      }
	    }
	    backgroundError[itype][ipt] = sqrt(sum / (nbins - 1));
	  
	    // Fill histograms with removed background
	    for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	      _histYieldBkgRemoved[itype][ipt]->fillBin(ibin, 
	        _histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt]);
	    }
	  
	    // Integrate near-side yield
	    unsigned int lowerBin = _histYield[itype][ipt]->binIndexAt(-0.7 + 0.02);
	    unsigned int upperBin = _histYield[itype][ipt]->binIndexAt( 0.7 - 0.02) + 1;
	    nbins = upperBin - lowerBin;
	    numberOfBins[itype][ipt][0] = nbins;
	    nearSide[itype][ipt] = _histYield[itype][ipt]->integralRange(lowerBin,
	      upperBin) - nbins * background[itype][ipt];
	    numberOfEntries[itype][ipt][0] = _histYield[itype][ipt]->integralRange(
	      lowerBin, upperBin) * _histTriggerCounter->bin(itype).sumW();
	  
	    // Integrate away-side yield
	    lowerBin = _histYield[itype][ipt]->binIndexAt(M_PI - 0.7 + 0.02);
	    upperBin = _histYield[itype][ipt]->binIndexAt(M_PI + 0.7 - 0.02) + 1;
	    nbins = upperBin - lowerBin;
	    numberOfBins[itype][ipt][1] = nbins;
	    awaySide[itype][ipt] = _histYield[itype][ipt]->integralRange(
	      lowerBin, upperBin) - nbins * background[itype][ipt];
	    numberOfEntries[itype][ipt][1] = 
	      _histYield[itype][ipt]->integralRange(lowerBin, upperBin) *
	      _histTriggerCounter->bin(itype).sumW();
	  
	  }
	}
      
	// Variables for IAA/ICP plots
	double dI = 0.0;
	int near = 0;
	int away = 1;
	double xval[PT_BINS] = { 3.5, 5.0, 7.0, 9.0 };
	double xerr[PT_BINS] = { 0.5, 1.0, 1.0, 1.0 };
      
	int types1[3] = {1, 2, 1};
	int types2[3] = {0, 0, 2};

	// Fill IAA/ICP plots for near side peak
	for (int ihist = 0; ihist < 3; ihist++) {
	  int type1 = types1[ihist];
	  int type2 = types2[ihist];
	  for (int ipt = 0; ipt < PT_BINS; ipt++) {
	    dI = scalingFactor[type1] * scalingFactor[type1] * numberOfEntries[type1][ipt][near] + 
	      scalingFactor[type2] * scalingFactor[type2] * numberOfEntries[type2][ipt][near] * 
	      nearSide[type1][ipt] * nearSide[type1][ipt] / (nearSide[type2][ipt] * nearSide[type2][ipt]) + 
	      numberOfBins[type1][ipt][near] * numberOfBins[type1][ipt][near] * backgroundError[type1][ipt] *
	      backgroundError[type1][ipt] + numberOfBins[type2][ipt][near] * numberOfBins[type2][ipt][near] * 
	      backgroundError[type2][ipt] * backgroundError[type2][ipt] * nearSide[type1][ipt] * 
	      nearSide[type1][ipt] / (nearSide[type2][ipt] * nearSide[type2][ipt]);

	    dI = sqrt(dI)/nearSide[type2][ipt];
	    _histIAA[ihist]->addPoint(xval[ipt], nearSide[type1][ipt] / nearSide[type2][ipt], xerr[ipt], dI);
	  }
	}
      
	// Fill IAA/ICP plots for away side peak
	for (int ihist = 0; ihist < 3; ihist++) {
	  int type1 = types1[ihist];
	  int type2 = types2[ihist];
	  for (int ipt = 0; ipt < PT_BINS; ipt++) {
	    dI = scalingFactor[type1] * scalingFactor[type1] * numberOfEntries[type1][ipt][away] + 
	      scalingFactor[type2] * scalingFactor[type2] * numberOfEntries[type2][ipt][away] * 
	      awaySide[type1][ipt] * awaySide[type1][ipt] / (awaySide[type2][ipt] * awaySide[type2][ipt]) + 
	      numberOfBins[type1][ipt][away] * numberOfBins[type1][ipt][away] * backgroundError[type1][ipt] *
	      backgroundError[type1][ipt] + numberOfBins[type2][ipt][away] * numberOfBins[type2][ipt][away] *
	      backgroundError[type2][ipt] * backgroundError[type2][ipt] * awaySide[type1][ipt] * 
	      awaySide[type1][ipt] / (awaySide[type2][ipt] * awaySide[type2][ipt]);

	    dI = sqrt(dI)/awaySide[type2][ipt];
	    _histIAA[ihist + 3]->addPoint(xval[ipt], awaySide[type1][ipt] / awaySide[type2][ipt], xerr[ipt], dI);
	  }
	}
      }
      
    }
    
  private:
    
    static const int PT_BINS = 4;
    static const int EVENT_TYPES = 3;

    Histo1DPtr _histYield[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histYieldBkgRemoved[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histTriggerCounter;
    Scatter2DPtr _histIAA[6];
    double pt_limits[5];
    int event_type;
    string event_string[EVENT_TYPES + 1];
    bool reentrant_flag;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I930312);
}
