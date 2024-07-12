// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

namespace Rivet {


  class CMS_2013_I1224539 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CMS_2013_I1224539()
      : Analysis("CMS_2013_I1224539"),
        _filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), fastjet::SelectorNHardest(3))),
        _trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.03))),
        _pruner(fastjet::Pruner(fastjet::cambridge_algorithm, 0.1, 0.5))
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Get options 
      WJET = true;
      ZJET = true;
      DIJET = true;
      if ( getOption("JMODE") == "W" ) {
	ZJET = false;
	DIJET = false;
      }
      if ( getOption("JMODE") == "Z" ) {
	WJET = false;
	DIJET = false;
      }
      if ( getOption("JMODE") == "DIJET" ) {
	WJET = false;
	ZJET = false;
      }


      FinalState fs(Cuts::abseta < 2.4);
      declare(fs, "FS");

      if (WJET) {
	// filling W+jet histos

	// Find W's with pT > 120, MET > 50
	WFinder wfinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 80*GeV, PID::ELECTRON, 50*GeV, 1000*GeV, 50.0*GeV,
			0.2, WFinder::ChargedLeptons::PROMPT, WFinder::ClusterPhotons::NODECAY, WFinder::AddPhotons::NO, WFinder::MassWindow::MT);
	declare(wfinder, "WFinder");

	// W+jet jet collections
	declare(FastJets(wfinder.remainingFinalState(), FastJets::ANTIKT, 0.7), "JetsAK7_wj");
	declare(FastJets(wfinder.remainingFinalState(), FastJets::CAM, 0.8), "JetsCA8_wj");
	declare(FastJets(wfinder.remainingFinalState(), FastJets::CAM, 1.2), "JetsCA12_wj");

	// Histograms
	/// @note These are 2D histos rendered into slices
	const int wjetsOffset = 51;
	for (size_t i = 0; i < N_PT_BINS_vj; ++i) {
	  book(_h_ungroomedJetMass_AK7_wj[i] ,wjetsOffset+i+1+0*N_PT_BINS_vj, 1, 1);
	  book(_h_filteredJetMass_AK7_wj[i] ,wjetsOffset+i+1+1*N_PT_BINS_vj, 1, 1);
	  book(_h_trimmedJetMass_AK7_wj[i] ,wjetsOffset+i+1+2*N_PT_BINS_vj, 1, 1);
	  book(_h_prunedJetMass_AK7_wj[i] ,wjetsOffset+i+1+3*N_PT_BINS_vj, 1, 1);
	  book(_h_prunedJetMass_CA8_wj[i] ,wjetsOffset+i+1+4*N_PT_BINS_vj, 1, 1);
	  if (i > 0) book(_h_filteredJetMass_CA12_wj[i] ,wjetsOffset+i+5*N_PT_BINS_vj, 1, 1);
	}
      }

      if (ZJET) {
	// filling Z+jet histos

	// Find Zs with pT > 120 GeV
	ZFinder zfinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 30*GeV, PID::ELECTRON, 80*GeV, 100*GeV,
			0.2, ZFinder::ClusterPhotons::NODECAY, ZFinder::AddPhotons::YES);
	declare(zfinder, "ZFinder");

	// Z+jet jet collections
	declare(FastJets(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.7), "JetsAK7_zj");
	declare(FastJets(zfinder.remainingFinalState(), FastJets::CAM, 0.8), "JetsCA8_zj");
	declare(FastJets(zfinder.remainingFinalState(), FastJets::CAM, 1.2), "JetsCA12_zj");

	// Histograms
	/// @note These are 2D histos rendered into slices
	const int zjetsOffset = 28;
	for (size_t i = 0; i < N_PT_BINS_vj; ++i ) {
	  book(_h_ungroomedJetMass_AK7_zj[i] ,zjetsOffset+i+1+0*N_PT_BINS_vj, 1, 1);
	  book(_h_filteredJetMass_AK7_zj[i] ,zjetsOffset+i+1+1*N_PT_BINS_vj,1,1);
	  book(_h_trimmedJetMass_AK7_zj[i] ,zjetsOffset+i+1+2*N_PT_BINS_vj,1,1);
	  book(_h_prunedJetMass_AK7_zj[i] ,zjetsOffset+i+1+3*N_PT_BINS_vj,1,1);
	  book(_h_prunedJetMass_CA8_zj[i] ,zjetsOffset+i+1+4*N_PT_BINS_vj,1,1);
	  if (i > 0) book(_h_filteredJetMass_CA12_zj[i] ,zjetsOffset+i+5*N_PT_BINS_vj,1,1);
	}

      }

      if (DIJET){

	// Jet collections
	declare(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");
	declare(FastJets(fs, FastJets::CAM, 0.8), "JetsCA8");
	declare(FastJets(fs, FastJets::CAM, 1.2), "JetsCA12");

	// Histograms
	for (size_t i = 0; i < N_PT_BINS_dj; ++i ) {
	  book(_h_ungroomedAvgJetMass_dj[i] ,i+1+0*N_PT_BINS_dj, 1, 1);
	  book(_h_filteredAvgJetMass_dj[i] ,i+1+1*N_PT_BINS_dj, 1, 1);
	  book(_h_trimmedAvgJetMass_dj[i] ,i+1+2*N_PT_BINS_dj, 1, 1);
	  book(_h_prunedAvgJetMass_dj[i] ,i+1+3*N_PT_BINS_dj, 1, 1);
	}

      }
    }

    bool isBackToBack_wj(const WFinder& wf, const fastjet::PseudoJet& psjet) {
      const FourMomentum w = wf.bosons()[0];
      const FourMomentum l1 = wf.constituentLeptons()[0];
      const FourMomentum l2 = wf.constituentNeutrinos()[0];
      /// @todo We should make FourMomentum know how to construct itself from a PseudoJet
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(w, jmom) > 2.0 && deltaR(l1, jmom) > 1.0 && deltaPhi(l2, jmom) > 0.4);
    }
    
    bool isBackToBack_zj(const ZFinder& zf, const fastjet::PseudoJet& psjet) {
      const FourMomentum& z = zf.bosons()[0].momentum();
      const FourMomentum& l1 = zf.constituents()[0].momentum();
      const FourMomentum& l2 = zf.constituents()[1].momentum();
      /// @todo We should make FourMomentum know how to construct itself from a PseudoJet
      const FourMomentum jmom(psjet.e(), psjet.px(), psjet.py(), psjet.pz());
      return (deltaPhi(z, jmom) > 2.0 && deltaR(l1, jmom) > 1.0 && deltaR(l2, jmom) > 1.0);
    }

    
    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin_vj(double ptJ) {
      const double ptBins_vj[N_PT_BINS_vj+1] = { 125.0, 150.0, 220.0, 300.0, 450.0 };
      for (size_t ibin = 0; ibin < N_PT_BINS_vj; ++ibin) {
	if (inRange(ptJ, ptBins_vj[ibin], ptBins_vj[ibin+1])) return ibin;
      }
      return N_PT_BINS_vj;
    }

    // Find the pT histogram bin index for value pt (in GeV), to hack a 2D histogram equivalent
    /// @todo Use a YODA axis/finder alg when available
    size_t findPtBin_jj(double ptJ) {
      const double ptBins_dj[N_PT_BINS_dj+1] = { 220.0, 300.0, 450.0, 500.0, 600.0, 800.0, 1000.0, 1500.0};
      for (size_t ibin = 0; ibin < N_PT_BINS_dj; ++ibin) {
        if (inRange(ptJ, ptBins_dj[ibin], ptBins_dj[ibin+1])) return ibin;
      }
      return N_PT_BINS_dj;
    }

  
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = 1.0;

      if (WJET) {

	// Get the W
	const WFinder& wfinder = apply<WFinder>(event, "WFinder");
	if (wfinder.bosons().size() == 1) {
	  const Particle w = wfinder.bosons()[0];
	  const Particle l = wfinder.constituentLeptons()[0];
	
	  // Require a fairly high-pT W and charged lepton
	  if (l.pT() >= 80*GeV && w.pT() >= 120*GeV) {
	
	      // Get the pseudojets.
	      const PseudoJets psjetsCA8_wj = apply<FastJets>(event, "JetsCA8_wj").pseudoJetsByPt( 50.0*GeV );
	      const PseudoJets psjetsCA12_wj = apply<FastJets>(event, "JetsCA12_wj").pseudoJetsByPt( 50.0*GeV );
	      
	      // AK7 jets
	      const PseudoJets psjetsAK7_wj = apply<FastJets>(event, "JetsAK7_wj").pseudoJetsByPt( 50.0*GeV );
	      if (!psjetsAK7_wj.empty()) {
		// Get the leading jet and make sure it's back-to-back with the W
		const fastjet::PseudoJet& j0 = psjetsAK7_wj[0];
		if (isBackToBack_wj(wfinder, j0)) {
		  const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
		  if (njetBin < N_PT_BINS_vj) {
		    fastjet::PseudoJet filtered0 = _filter(j0);
		    fastjet::PseudoJet trimmed0 = _trimmer(j0);
		    fastjet::PseudoJet pruned0 = _pruner(j0);
		    _h_ungroomedJetMass_AK7_wj[njetBin]->fill(j0.m()/GeV, weight);
		    _h_filteredJetMass_AK7_wj[njetBin]->fill(filtered0.m()/GeV, weight);
		    _h_trimmedJetMass_AK7_wj[njetBin]->fill(trimmed0.m()/GeV, weight);
		    _h_prunedJetMass_AK7_wj[njetBin]->fill(pruned0.m()/GeV, weight);
		  }
		}
	      }
	      
	      // CA8 jets
	      if (!psjetsCA8_wj.empty()) {
		// Get the leading jet and make sure it's back-to-back with the W
		const fastjet::PseudoJet& j0 = psjetsCA8_wj[0];
		if (isBackToBack_wj(wfinder, j0)) {
		  const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
		  if (njetBin < N_PT_BINS_vj) {
		    fastjet::PseudoJet pruned0 = _pruner(j0);
		    _h_prunedJetMass_CA8_wj[njetBin]->fill(pruned0.m()/GeV, weight);
		  }
		}
	      }
	      
	      // CA12 jets
	      if (!psjetsCA12_wj.empty()) {
		// Get the leading jet and make sure it's back-to-back with the W
		const fastjet::PseudoJet& j0 = psjetsCA12_wj[0];
		if (isBackToBack_wj(wfinder, j0)) {
		  const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
		  if (njetBin < N_PT_BINS_vj&&njetBin>0) {
		    fastjet::PseudoJet filtered0 = _filter(j0);
		    _h_filteredJetMass_CA12_wj[njetBin]->fill( filtered0.m() / GeV, weight);
		  }
		}
	      }
	    }
	}
      }

      if (ZJET) {
	// Get the Z
	const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
	if (zfinder.bosons().size() == 1) {
	  const Particle& z = zfinder.bosons()[0];
	  if (z.constituents().size() < 2) {
	    MSG_WARNING("Found a Z with less than 2 constituents.");
	    vetoEvent;
	  }
	  const Particle& l1 = z.constituents()[0];
	  const Particle& l2 = z.constituents()[1];
	  MSG_DEBUG(l1.pT() << " " << l2.pT());
	  assert(&l1 != &l2);
	  
	  // Require a high-pT Z (and constituents)
	  if (l1.pT() >= 30*GeV && l2.pT() >= 30*GeV && z.pT() >= 120*GeV) {

	    // AK7 jets
	    const PseudoJets& psjetsAK7_zj = apply<FastJets>(event, "JetsAK7_zj").pseudoJetsByPt(50.0*GeV);
	    if (!psjetsAK7_zj.empty()) {
	      // Get the leading jet and make sure it's back-to-back with the Z
	      const fastjet::PseudoJet& j0 = psjetsAK7_zj[0];
	      if (isBackToBack_zj(zfinder, j0)) {
		const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
		if (njetBin < N_PT_BINS_vj) {
		  fastjet::PseudoJet filtered0 = _filter(j0);
		  fastjet::PseudoJet trimmed0 = _trimmer(j0);
		  fastjet::PseudoJet pruned0 = _pruner(j0);
		  _h_ungroomedJetMass_AK7_zj[njetBin]->fill(j0.m()/GeV, weight);
		  _h_filteredJetMass_AK7_zj[njetBin]->fill(filtered0.m()/GeV, weight);
		  _h_trimmedJetMass_AK7_zj[njetBin]->fill(trimmed0.m()/GeV, weight);
		  _h_prunedJetMass_AK7_zj[njetBin]->fill(pruned0.m()/GeV, weight);
		}
	      }
	    }
	    
	    // CA8 jets
	    const PseudoJets& psjetsCA8_zj = apply<FastJets>(event, "JetsCA8_zj").pseudoJetsByPt(50.0*GeV);
	    if (!psjetsCA8_zj.empty()) {
	      // Get the leading jet and make sure it's back-to-back with the Z
	      const fastjet::PseudoJet& j0 = psjetsCA8_zj[0];
	      if (isBackToBack_zj(zfinder, j0)) {
		const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
		if (njetBin < N_PT_BINS_vj) {
		  fastjet::PseudoJet pruned0 = _pruner(j0);
		  _h_prunedJetMass_CA8_zj[njetBin]->fill(pruned0.m()/GeV, weight);
		}
	      }
	    }
	    
	    // CA12 jets
	    const PseudoJets& psjetsCA12_zj = apply<FastJets>(event, "JetsCA12_zj").pseudoJetsByPt(50.0*GeV);
	    if (!psjetsCA12_zj.empty()) {
	      // Get the leading jet and make sure it's back-to-back with the Z
	      const fastjet::PseudoJet& j0 = psjetsCA12_zj[0];
	      if (isBackToBack_zj(zfinder, j0)) {
		const size_t njetBin = findPtBin_vj(j0.pt()/GeV);
		if (njetBin>0 && njetBin < N_PT_BINS_vj) {
		  fastjet::PseudoJet filtered0 = _filter(j0);
		  _h_filteredJetMass_CA12_zj[njetBin]->fill( filtered0.m() / GeV, weight);
		}
	      }
	    }
	  }
	}
      }

      if (DIJET) {
	// Look at events with >= 2 jets
	const PseudoJets& psjetsAK7 = apply<FastJets>(event, "JetsAK7").pseudoJetsByPt( 50.0*GeV );
	if (psjetsAK7.size() >= 2) {
	
	  // Get the leading two jets and find their average pT
	  const fastjet::PseudoJet& j0 = psjetsAK7[0];
	  const fastjet::PseudoJet& j1 = psjetsAK7[1];
	  double ptAvg = 0.5 * (j0.pt() + j1.pt());
	  
	  // Find the appropriate mean pT bin and escape if needed
	  const size_t njetBin = findPtBin_jj(ptAvg/GeV);
	  if (njetBin < N_PT_BINS_dj) {
	  
	    // Now run the substructure algs...
	    fastjet::PseudoJet filtered0 = _filter(j0);
	    fastjet::PseudoJet filtered1 = _filter(j1);
	    fastjet::PseudoJet trimmed0 = _trimmer(j0);
	    fastjet::PseudoJet trimmed1 = _trimmer(j1);
	    fastjet::PseudoJet pruned0 = _pruner(j0);
	    fastjet::PseudoJet pruned1 = _pruner(j1);
	    
	    // ... and fill the histograms
	    _h_ungroomedAvgJetMass_dj[njetBin]->fill(0.5*(j0.m() + j1.m())/GeV, weight);
	    _h_filteredAvgJetMass_dj[njetBin]->fill(0.5*(filtered0.m() + filtered1.m())/GeV, weight);
	    _h_trimmedAvgJetMass_dj[njetBin]->fill(0.5*(trimmed0.m() + trimmed1.m())/GeV, weight);
	    _h_prunedAvgJetMass_dj[njetBin]->fill(0.5*(pruned0.m() + pruned1.m())/GeV, weight);
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double normalizationVal = 1000;


      for (size_t i = 0; i < N_PT_BINS_vj; ++i) {
	if (WJET) {
	  normalize(_h_ungroomedJetMass_AK7_wj[i], normalizationVal);
	  normalize(_h_filteredJetMass_AK7_wj[i], normalizationVal);
	  normalize(_h_trimmedJetMass_AK7_wj[i], normalizationVal);
	  normalize(_h_prunedJetMass_AK7_wj[i], normalizationVal);
	  normalize(_h_prunedJetMass_CA8_wj[i], normalizationVal);
	  if (i > 0) normalize( _h_filteredJetMass_CA12_wj[i], normalizationVal);
	}
	if (ZJET) {
	  normalize( _h_ungroomedJetMass_AK7_zj[i], normalizationVal);
	  normalize( _h_filteredJetMass_AK7_zj[i], normalizationVal);
	  normalize( _h_trimmedJetMass_AK7_zj[i], normalizationVal);
	  normalize( _h_prunedJetMass_AK7_zj[i], normalizationVal);
	  normalize( _h_prunedJetMass_CA8_zj[i], normalizationVal);
	  if (i > 0) normalize( _h_filteredJetMass_CA12_zj[i], normalizationVal);
	}
      }
      if (DIJET) {
	for (size_t i = 0; i < N_PT_BINS_dj; ++i) {
	  normalize(_h_ungroomedAvgJetMass_dj[i], normalizationVal);
	  normalize(_h_filteredAvgJetMass_dj[i], normalizationVal);
	  normalize(_h_trimmedAvgJetMass_dj[i], normalizationVal);
	  normalize(_h_prunedAvgJetMass_dj[i], normalizationVal);
	}
      }
    }

    //@}

  protected:

    bool WJET, ZJET, DIJET;

  private:

    /// @name FastJet grooming tools (configured in constructor init list)
    //@{
    const fastjet::Filter _filter;
    const fastjet::Filter _trimmer;
    const fastjet::Pruner _pruner;
    //@}


    /// @name Histograms
    //@{
    enum BINS_vj { PT_125_150_vj=0, PT_150_220_vj, PT_220_300_vj, PT_300_450_vj, N_PT_BINS_vj };
    // W+jet
    Histo1DPtr _h_ungroomedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_trimmedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_AK7_wj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_CA8_wj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_CA12_wj[N_PT_BINS_vj];
    //Z+jet
    Histo1DPtr _h_ungroomedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_trimmedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_AK7_zj[N_PT_BINS_vj];
    Histo1DPtr _h_prunedJetMass_CA8_zj[N_PT_BINS_vj];
    Histo1DPtr _h_filteredJetMass_CA12_zj[N_PT_BINS_vj];
    // DIJET
    enum BINS_dj { PT_220_300_dj=0, PT_300_450_dj, PT_450_500_dj, PT_500_600_dj,
                   PT_600_800_dj, PT_800_1000_dj, PT_1000_1500_dj, N_PT_BINS_dj };
    Histo1DPtr _h_ungroomedJet0pt, _h_ungroomedJet1pt;
    Histo1DPtr _h_ungroomedAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_filteredAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_trimmedAvgJetMass_dj[N_PT_BINS_dj];
    Histo1DPtr _h_prunedAvgJetMass_dj[N_PT_BINS_dj];
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CMS_2013_I1224539);

}
