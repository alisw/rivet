// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief multiplicities in u, d, s events
  class OPAL_2001_I536266 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2001_I536266);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");
      book(_cDown   , "/TMP/CDOWN"   );
      book(_cUp     , "/TMP/CUP"     );
      book(_cStrange, "/TMP/CSTRANGE");
      book(_wDown   , "/TMP/WDOWN"   );
      book(_wUp     , "/TMP/WUP"     );
      book(_wStrange, "/TMP/WSTRANGE");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;


      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      const size_t numParticles = cfs.particles().size();
      switch (flavour) {
      case 1:
	_wDown->fill();
        _cDown->fill(numParticles);
        break;
      case 2:
        _wUp->fill();
	_cUp->fill(numParticles);
        break;
      case 3:
        _wStrange->fill();
        _cStrange->fill(numParticles);
        break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // calculate the averages and ratios
      if(_wUp      ->effNumEntries()!=0. ) scale( _cUp     , 1./ *_wUp);
      if(_wDown    ->effNumEntries()!=0. ) scale( _cDown   , 1./ *_wDown);
      if(_wStrange ->effNumEntries()!=0. ) scale( _cStrange, 1./ *_wStrange);
      for(unsigned int ix=1;ix<3;++ix) {
	for(unsigned int iy=1;iy<4;++iy) {
	  double val;
	  std::pair<double,double> errs;
	  if(ix==1) {
	    CounterPtr cTemp;
	    if(iy==1)      cTemp = _cUp;
	    else if(iy==2) cTemp = _cDown;
	    else if(iy==3) cTemp = _cStrange;
	    val  = cTemp->val();
	    errs = make_pair(cTemp->err(),cTemp->err()); 
	  }
	  else {
	    Scatter1D temp;
	    if(iy==1)      temp = *_cUp     / *_cDown;
	    else if(iy==2) temp = *_cStrange/ *_cDown;
	    else if(iy==3) temp = *_cStrange/ *_cUp  ;
	    val  = temp.points()[0].x();
	    errs = temp.points()[0].xErrs(); 
	  }
	  Scatter2DPtr  mult;
	  book(mult, ix, 1, iy);
	  mult->addPoint(45.6, val, make_pair(0.5,0.5), errs);
	}
      }
    }

    //@}

    /// @name Multiplicities
    //@{
    CounterPtr _cDown;
    CounterPtr _cUp;
    CounterPtr _cStrange;
    //@}

    /// @name Weights
    //@{
    CounterPtr _wDown;
    CounterPtr _wUp;
    CounterPtr _wStrange;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OPAL_2001_I536266);


}
