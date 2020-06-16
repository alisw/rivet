// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D_s D_s* spectrum
  class BABAR_2002_I582184 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2002_I582184);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      // rates
      book(_c_Ds_on     ,2,1,1);
      book(_c_Ds_off    ,1,1,1);
      book(_c_DsStar_on ,2,1,2);
      book(_c_DsStar_off,1,1,2);
      book(_w_ups ,"/TMP/w_ups" );
      // dists
      book(_h_Ds_on     ,5,1,1);
      book(_h_Ds_off    ,3,1,2);
      book(_h_DsStar_on ,5,1,2);
      book(_h_DsStar_off,4,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilon(4S) among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      bool cont = ufs.particles(Cuts::pid==300553).empty();
      if(!cont) _w_ups ->fill();
      for(const Particle &p : ufs.particles(Cuts::abspid==431 or Cuts::abspid==433)) {
	double mom=p.momentum().p3().mod();
	if(cont) {
	  if(p.abspid()==431) {
	    _h_Ds_off->fill(mom);
	    _c_Ds_off->fill(0.5);
	  }
	  else {
	    _h_DsStar_off->fill(mom);
	    _c_DsStar_off->fill(0.5);
	  }
	}
	else {
	  if(p.abspid()==431) {
	    _h_Ds_on->fill(mom);
	    _c_Ds_on->fill(0.5);
	  }
	  else {
	    _h_DsStar_on->fill(mom);
	    _c_DsStar_on->fill(0.5);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_Ds_on     );
      normalize(_h_Ds_off    );
      normalize(_h_DsStar_on );
      normalize(_h_DsStar_off);
      if(_w_ups->val()!=0) {
	scale(_c_Ds_on    ,0.5/ *_w_ups);
	scale(_c_DsStar_on,0.5/ *_w_ups);
      }
      scale(_c_Ds_off    , crossSection()/sumOfWeights()/picobarn);
      scale(_c_DsStar_off, crossSection()/sumOfWeights()/picobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Ds_on,_h_Ds_off,_h_DsStar_on,_h_DsStar_off;
    Histo1DPtr _c_Ds_on,_c_Ds_off,_c_DsStar_on,_c_DsStar_off;
    CounterPtr _w_ups;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2002_I582184);

}
