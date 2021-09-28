// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, K and p spectra at 34 and 44 GeV
  class TASSO_1989_I267755 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1989_I267755);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      // Book histograms
      _iHist=-1;
      if(fuzzyEquals(sqrtS()/GeV, 34., 1e-3)) {
	_iHist = 0;
      }
      else if (fuzzyEquals(sqrtS()/GeV, 44., 1e-3)) {
	_iHist = 1;
      }
      else
	MSG_ERROR("Beam energy not supported!");

      book(_h_x_pi, 3*_iHist+7,1,1);
      book(_h_x_K , 3*_iHist+8,1,1);
      book(_h_x_p , 3*_iHist+9,1,1);
      if(_iHist==1) book(_h_x_pi0,13,1,1);
      book(_n_pi,"TMP/n_pi",refData(3*_iHist+1,1,1));
      book(_d_pi,"TMP/d_pi",refData(3*_iHist+1,1,1));
      book(_n_K ,"TMP/n_K" ,refData(3*_iHist+2,1,1));
      book(_d_K ,"TMP/d_K" ,refData(3*_iHist+2,1,1));
      book(_n_p ,"TMP/n_p" ,refData(3*_iHist+3,1,1));
      book(_d_p ,"TMP/d_p" ,refData(3*_iHist+3,1,1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      
      for (const Particle& p : fs.particles()) {
      	double xP = p.p3().mod()/meanBeamMom;
      	_d_pi->fill(xP);
      	_d_K ->fill(xP);
      	_d_p ->fill(xP);
      	if(abs(p.pid())==211) {
      	  _h_x_pi->fill(xP);
      	  _n_pi  ->fill(xP);
      	}
      	else if(abs(p.pid())==321) {
      	  _h_x_K->fill(xP);
      	  _n_K  ->fill(xP);
      	}
      	else if(abs(p.pid())==2212) {
      	  _h_x_p->fill(xP);
      	  _n_p  ->fill(xP);
      	}
      }
      if(_h_x_pi0) {
      	for(const Particle & p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==111)) {
      	  double xP = p.p3().mod()/meanBeamMom;
      	  _h_x_pi0->fill(xP);
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_x_pi , 1./sumOfWeights());
      scale(_h_x_K  , 1./sumOfWeights());
      scale(_h_x_p  , 1./sumOfWeights());
      if(_h_x_pi0) scale(_h_x_pi0, 1./sumOfWeights());
      Scatter2DPtr temp1,temp2,temp3;
      book(temp1,3*_iHist+1,1,1);
      book(temp2,3*_iHist+2,1,1);
      book(temp3,3*_iHist+3,1,1);
      
      divide(_n_pi,_d_pi, temp1);
      divide(_n_K ,_d_K , temp2);
      divide(_n_p ,_d_p , temp3);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_p_pi, _h_x_pi, _h_p_K, _h_x_K, _h_p_p, _h_x_p, _h_x_pi0;
    Histo1DPtr _n_pi,_d_pi,_n_K,_d_K,_n_p,_d_p;
    int _iHist;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1989_I267755);


}
