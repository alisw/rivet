// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Charged particle spectra at 5.2, 6.5 and 29 GeV
  class MARKII_1982_I178416 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1982_I178416);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState fs;
      declare(fs, "FS");
      unsigned int iloc(0);
      if(isCompatibleWithSqrtS(5.2))
	iloc = 1;
      else if(isCompatibleWithSqrtS(6.5))
	iloc = 2;
      else if(isCompatibleWithSqrtS(29.0))
	iloc = 3;
      else
        MSG_ERROR("Beam energy incompatible with analysis.");
      assert(iloc!=0);
      book(_h_x,1,1,iloc);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      if(fs.particles().size()==2 &&
	 abs(fs.particles()[0].pid())==13 &&
	 abs(fs.particles()[1].pid())==13) vetoEvent;
      for (const Particle& p : fs.particles()) {
	const Vector3 mom3 = p.p3();
	double pp = mom3.mod();
	double x = 2.*pp/sqrtS();
	_h_x->fill(x);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_x,crossSection()*sqr(sqrtS())/sumOfWeights()/microbarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1982_I178416);

}
