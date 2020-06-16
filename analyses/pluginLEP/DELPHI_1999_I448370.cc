// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class DELPHI_1999_I448370 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_1999_I448370);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");
      
      // Book histograms
       book(_h_F_T       , 1, 1, 1);
       book(_h_F_L       , 2, 1, 1);
       book(_h_F_A       , 3, 1, 1);
       book(_h_F_TL      , 4, 1, 1);
       book(_h_F_T_total , 5, 1, 1);
       book(_h_F_L_total , 5, 1, 2);
       book(_h_F_TL_total, 5, 1, 3);
       book(_h_b_F_T     , 7, 1, 1);
       book(_h_b_F_L     , 7, 1, 2);
       book(_h_light_F_T , 8, 1, 1);
       book(_h_light_F_L , 8, 1, 2);
       book(_n_bottom    , 9, 1, 1);
       book(_n_light     , 9, 1, 2);
      
       book(_c_light , "/TMP/wLight");
       book(_c_bottom, "/TMP/wBottom");
       book(_c_total , "/TMP/wTotal");

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
      if     (flavour==5) _c_bottom->fill();
      else if(flavour!=4) _c_light ->fill();
      _c_total->fill();
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();

      // loop over charged particles
      double v=0.8, v2=sqr(v), v5=v2*v2*v;
      for (const Particle& p : fs.particles()) {
        double xp = p.p3().mod()/meanBeamMom;
	double ctheta = axis.dot(p.momentum().p3().unit());
	if(abs(ctheta)<0.8) {
	  double WT = 0.5 /v5*(5.*sqr(ctheta)*(3.-v2)-v2*(5.-3.*v2));
	  double WL = 0.25/v5*(v2*(5.+3.*v2)-5.*sqr(ctheta)*(3.+v2));
	  double WA = 2.*ctheta/v2/v;
	  _h_F_T  ->fill(xp,WT);
	  _h_F_L  ->fill(xp,WL);
	  _h_F_TL ->fill(xp,(WL+WT));
	  _h_F_T_total ->fill(xp,0.5*xp*WT);
	  _h_F_L_total ->fill(xp,0.5*xp*WL);
	  _h_F_TL_total->fill(xp,0.5*xp*(WL+WT));
	  if(p.charge()>0) {
	    _h_F_A->fill(xp, WA);
	  }
	  else {
	    _h_F_A->fill(xp,-WA);
	  }
	  if(flavour==5) {
	    _h_b_F_T  ->fill(xp,WT);
	    _h_b_F_L  ->fill(xp,WL);
	  }
	  else if(flavour!=4) {
	    _h_light_F_T->fill(xp,WT);
	    _h_light_F_L->fill(xp,WL);
	  }
	}
      }

      if(flavour==5)
	_n_bottom->fill(91.2,fs.particles().size());
      else if(flavour!=4)
	_n_light->fill(91.2,fs.particles().size());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_F_T ,1./ *_c_total);
      scale(_h_F_L ,1./ *_c_total);
      scale(_h_F_A ,1./ *_c_total);
      scale(_h_F_TL,1./ *_c_total);
      {Scatter2DPtr temp; divide(_h_F_L, _h_F_T , book(temp, 6, 1, 1));}
      {Scatter2DPtr temp; divide(_h_F_L, _h_F_TL, book(temp, 6, 1, 2));}
      scale(_h_F_T_total ,1./ *_c_total);
      scale(_h_F_L_total ,1./ *_c_total);
      scale(_h_F_TL_total,1./ *_c_total);
      scale(_h_b_F_T ,1./ *_c_bottom);
      scale(_h_b_F_L ,1./ *_c_bottom);
      scale(_h_light_F_T ,1./ *_c_light);
      scale(_h_light_F_L ,1./ *_c_light);
      scale(_n_bottom ,1./ *_c_bottom);
      scale(_n_light ,1./ *_c_light);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_F_T,_h_F_L,_h_F_A,_h_F_TL,_h_b_F_T,_h_b_F_L,_h_light_F_T,_h_light_F_L,_n_light,_n_bottom;
    Histo1DPtr _h_F_T_total,_h_F_L_total,_h_F_TL_total;
    CounterPtr _c_light,_c_bottom,_c_total;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_1999_I448370);


}
