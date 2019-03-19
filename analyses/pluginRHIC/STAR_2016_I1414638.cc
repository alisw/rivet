// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Tools/Percentile.hh"
#include <complex>

namespace Rivet {
  /// @brief Centrality projection for STAR AuAu.
  class CentralMultiplicityCentrality : public SingleValueProjection {
  
  public:
    
    /// Constructor
    CentralMultiplicityCentrality() {
      declare(ChargedFinalState(Cuts::abseta < 0.5), "FSCentralMultCent");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(CentralMultiplicityCentrality);
  
  protected:

    /// Perform the projection
    void project(const Event& e) {
      clear();
      double estimate = 
        apply<FinalState>(e, "FSCentralMultCent").particles().size();
      set(estimate);
    }

    /// Compare projections.
    int compare(const Projection& p) const {
      return mkNamedPCmp(p, "FSCentralMultCent");
    }
  
  };

  /// @brief Centrality calibration analysis skeleton. Should be put in a
  /// different file.
  class STAR_BES_CALIB : public Analysis {
  
  public:
    
    STAR_BES_CALIB () : Analysis("STAR_BES_CALIB") {};
  /// Book histograms and initialise projections before the run
  void init() {

    // One projection for the actual observable, and one for the
    // generated impact parameter.
    declare(CentralMultiplicityCentrality(), "Centrality");
    declare(ImpactParameterProjection(), "IMP");

    // The calibration histogram:
    _calib = bookHisto1D("CMULT", 100, 0.0, 2000.0);

    // If histogram was pre-loaded, the calibration is done.
    _done = ( _calib->numEntries() > 0 );

    // The alternative histogram based on impact parameter. Note that
    // it MUST be named the same as the histogram for the experimental
    // observable with an added _IMP suffix for the Pecentile<>
    // binning to work properly.
    _impcalib = bookHisto1D("CMULT_IMP", 400, 0.0, 20.0);


  }
  
  /// Perform the per-event analysis
  void analyze(const Event& event) {

    if ( _done ) return;
    
    const double weight = event.weight();

    // The alternative centrality based on generated impact
    // parameter, assumes that the generator does not describe the
    // full final state, and should therefore be filled even if the
    // event is not triggered.
    _impcalib->fill(apply<SingleValueProjection>(event, "IMP")(), weight);

    _calib->fill(apply<SingleValueProjection>(event, "Centrality")(), weight);

  }
  
  /// Finalize
  void finalize() {

    _calib->normalize();
    _impcalib->normalize();

  }

private:

  /// The calibration histograms.
  Histo1DPtr _calib;
  Histo1DPtr _impcalib;

  /// Safeguard from adding to a pre-loaded histogram.
  bool _done;
  
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_BES_CALIB);

  /// @brief Third harmonic of azimuthal correlations in Au+Au collisions 
  //  at different COM energies.
  class STAR_2016_I1414638 : public Analysis {
  public: 
    /// Constructor	  
    STAR_2016_I1414638 () : Analysis("STAR_2016_I1414638") {}

    void init() {
	/// Projections
	/// The centrality projection. 
	declareCentrality(CentralMultiplicityCentrality(), 
	  "STAR_BES_CALIB", "CMULT", "CMULT");

	/// The observed particles.
    	declare(ChargedFinalState(Cuts::abseta < 1.0 &&
			       	   Cuts::pT > 0.2*GeV), "CFS");
	/// The centrality bins
	centralityBins = {5., 10, 20, 30, 40, 50, 60, 70, 80};
	/// The corresponding histograms for the different analysis energies.
	vector<double> energies = {7.7, 11.5, 14.5, 19.6, 27.0, 39.0, 
	  62.4, 200.0};
	int energy = -1;
	for (int i = 0, N = energies.size(); i < N; ++i) {
	  if (fuzzyEquals(sqrtS()/197./GeV,energies[i],1E-1)) energy = i;
	}
	if (energy == -1) MSG_ERROR("Incompatible beam energy!");
	for (int i = 0; i < 9; ++i)
	  h_v32[centralityBins[i]] = bookProfile1D(1 + i + 9 * energy, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = applyProjection<ChargedFinalState>(event, "CFS");
      // Require at least two charged particles for the analysis to make sense.
      // No further triggers are described in the paper.
      const Particles& particles = cfs.particles();
      if (particles.size() < 2) return;
      // The centrality projection
      const CentralityProjection& cent = apply<CentralityProjection>(event,"CMULT");
      const double c = cent();
      // Find the correct histogram to fill.
      auto hItr = h_v32.upper_bound(c);
      if (hItr == h_v32.end()) return;
      for(int i = 0, N = particles.size(); i < N; ++i){
        for(int j = i + 1; j < N; ++j){
          const double eta1 = particles[i].eta();
	  const double eta2 = particles[j].eta();
	  if(eta1 * eta2 < 0){
            const double deltaPhi = abs(particles[i].phi() - particles[j].phi());
            // Fill profile with v_2(2)^2 from eq. (1) in the paper.
	    hItr->second->fill(abs(eta1 - eta2), cos(3.*deltaPhi),  weight);
          }
        }
       }
    }
    
    /// Normalise histograms etc., after the run
    void finalize() {
    
    }

    //@}


    /// @name Histograms
    //@{
    // The centralities.
    vector<double> centralityBins;
    // The histograms.
    map<double, Profile1DPtr> h_v32;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(STAR_2016_I1414638);


}
