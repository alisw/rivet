// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TPC_1988_I262143 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TPC_1988_I262143);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_h_z_pi ,1, 1, 1);
      book(_h_z_K  ,1, 1, 2);
      book(_h_z_p  ,1, 1, 3);
      book(_h_z_all,1, 1, 4);
      
      book(_h_z2_pi, 5, 1, 1);
      book(_h_z2_K , 5, 1, 2);
      book(_h_z2_p , 5, 1, 3);
      
      book(_n_pi,"TMP/n_pi", refData(6,1,1));
      book(_n_K ,"TMP/n_K" , refData(6,1,2));
      book(_n_p ,"TMP/n_p" , refData(6,1,3));
      book(_d_pi,"TMP/d_pi", refData(6,1,1));
      book(_d_K ,"TMP/d_K" , refData(6,1,2));
      book(_d_p ,"TMP/d_p" , refData(6,1,3));
      book(_n2_K,"TMP/n2_K", refData(7,1,1));
      book(_n2_p,"TMP/n2_p", refData(7,1,2));
      book(_n3_p,"TMP/n3_p", refData(7,1,3));
      book(_d2_K,"TMP/d2_K", refData(7,1,1));
      book(_d2_p,"TMP/d2_p", refData(7,1,2));
      book(_d3_p,"TMP/d3_p", refData(7,1,3));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
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
	_h_z_all->fill(xP);
	_d_pi->fill(xP);
	_d_K ->fill(xP);
	_d_p ->fill(xP);
	int id = abs(p.pid());
	if(id==211) {
	  _h_z_pi->fill(xP);
	  _h_z2_pi->fill(xP, xP);
	  _n_pi ->fill(xP, 100.);
	  _d2_K->fill(xP);
	  _d2_p->fill(xP);
	  _d3_p->fill(xP);
	}
	else if(id==321) {
	  _h_z_K ->fill(xP);
	  _h_z2_K ->fill(xP, xP);
	  _n_K ->fill(xP, 100.);
	  _n2_K->fill(xP);
	  _d3_p->fill(xP);
	}
	else if(id==2212) {
	  _h_z_p ->fill(xP);
	  _h_z2_p ->fill(xP, xP);
	  _n_p  ->fill(xP, 100.);
	  _n2_p->fill(xP);
	  _n3_p->fill(xP);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_z_all,1./sumOfWeights());
      scale(_h_z_pi ,1./sumOfWeights());
      scale(_h_z_K  ,1./sumOfWeights());
      scale(_h_z_p  ,1./sumOfWeights());
      scale(_h_z2_pi ,1./sumOfWeights());
      scale(_h_z2_K  ,1./sumOfWeights());
      scale(_h_z2_p  ,1./sumOfWeights());
      Scatter2DPtr temp1,temp2,temp3,temp4,temp5,temp6;
      book(temp1,6, 1, 1);
      book(temp2,6, 1, 2);
      book(temp3,6, 1, 3);
      book(temp4,7, 1, 1);
      book(temp5,7, 1, 2);
      book(temp6,7, 1, 3);
      divide(_n_pi,_d_pi, temp1);
      divide(_n_K ,_d_K , temp2);
      divide(_n_p ,_d_p , temp3);
      divide(_n2_K,_d2_K, temp4);
      divide(_n2_p,_d2_p, temp5);
      divide(_n3_p,_d3_p, temp6);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_z_all,_h_z_pi,_h_z_K,_h_z_p;
    Histo1DPtr _h_z2_pi,_h_z2_K,_h_z2_p;
    Histo1DPtr _n_pi,_n_K,_n_p,_d_pi,_d_K,_d_p;
    Histo1DPtr _n2_K,_n2_p,_n3_p,_d2_K,_d2_p,_d3_p;
    //@}2


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TPC_1988_I262143);


}
