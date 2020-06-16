// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Charged multiplicity for different numbers of final state jets
  class DELPHI_1992_I334948 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1992_I334948);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(FastJets(cfs, FastJets::JADE, 0.7), "Jets");

      // Book histograms
      for(unsigned int ih=0;ih<3;++ih) {
	for(unsigned int iy=0;iy<3;++iy) {
	  book(_h_mult[ih][iy],ih+1,1,iy+1);
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");
      const FastJets& jets = apply<FastJets>(event, "Jets");
      if (jets.clusterSeq()) {
	vector<double> ycut = {0.01,0.02,0.04};
	for (unsigned int ih=0;ih<3;++ih) {
	  int nbin = jets.clusterSeq()->n_exclusive_jets_ycut(ycut[ih])-2;
	  if(nbin<0 || nbin>2) continue;
	  _h_mult[ih][nbin]->fill(numParticles);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      for(unsigned int ih=0;ih<3;++ih) {
	for(unsigned int iy=0;iy<3;++iy) {
	  normalize(_h_mult[ih][iy],2000.);
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_mult[3][3];
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1992_I334948);


}
