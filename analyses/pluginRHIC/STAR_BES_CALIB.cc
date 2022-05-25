// -*- C++ -*-
#include <iostream>
#include <string>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ImpactParameterProjection.hh"
#include "Rivet/Projections/SingleValueProjection.hh"
#include "Rivet/Tools/Percentile.hh"
#include "Rivet/Tools/RHICCommon.hh"

namespace Rivet {

class STAR_BES_CALIB : public Analysis {
       public:
	STAR_BES_CALIB() : Analysis("STAR_BES_CALIB") {}
	void init() {
		declare(STAR_BES_Centrality(), "Centrality");
		declare(ImpactParameterProjection(), "IMP");
		// The calibrationhistogram:
		book(_calib, "CMULT", 100, 0.0, 200.0);

		// The alternative histogram based on impact parameter. Note
		// that it MUST be named the same as the histogram for the
		// experimental observable with an added _IMP suffix for the
		// Pecentile<> binning to work properly.
		book(_impcalib, "CMULT_IMP", 400, 0.0, 20.0);
	}

	void analyze(const Event& event) {

		_impcalib->fill(apply<SingleValueProjection>(event, "IMP")());

		_calib->fill(apply<STAR_BES_Centrality>(event, "Centrality")());
	}

	void finalize() {
		_calib->normalize();
		_impcalib->normalize();
	}

       private:
	/// The calibration histograms.
	Histo1DPtr _calib;
	Histo1DPtr _impcalib;

};

RIVET_DECLARE_PLUGIN(STAR_BES_CALIB);

}

