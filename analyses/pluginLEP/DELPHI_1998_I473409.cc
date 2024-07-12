// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief flavour seperate pi,K,p spectra
  class DELPHI_1998_I473409 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1998_I473409);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF"); 
      // Book histograms
      book(_h_all_pi, "TMP/h_all_pi",refData( 4,1,1));
      book(_h_all_K , "TMP/h_all_K ",refData( 5,1,1));
      book(_h_all_p , "TMP/h_all_p ",refData( 6,1,1));
      book(_h_all_Kp, "TMP/h_all_Kp",refData( 7,1,1));
      book(_d_all   , "TMP/d_all   ",refData( 4,1,1));
      			           
      book(_h_bot_pi, "TMP/h_bot_pi",refData( 8,1,1));
      book(_h_bot_K , "TMP/h_bot_K ",refData( 9,1,1));
      book(_h_bot_p , "TMP/h_bot_p ",refData(10,1,1));
      book(_h_bot_Kp, "TMP/h_bot_Kp",refData(11,1,1));
      book(_d_bot   , "TMP/d_bot   ",refData( 8,1,1));
      			           
      book(_h_lgt_pi, "TMP/h_lgt_pi",refData(12,1,1));
      book(_h_lgt_K , "TMP/h_lgt_K ",refData(13,1,1));
      book(_h_lgt_p , "TMP/h_lgt_p ",refData(14,1,1));
      book(_h_lgt_Kp, "TMP/h_lgt_Kp",refData(15,1,1));
      book(_d_lgt   , "TMP/d_lgt   ",refData(12,1,1));

      book(_h_all_ch_p, 16,1,1);
      book(_h_all_ch_x, 17,1,1);
      book(_h_all_pi_p, 18,1,1);
      book(_h_all_pi_x, 19,1,1);
      book(_h_all_K_p , 20,1,1);
      book(_h_all_k_x , 21,1,1);
      book(_h_all_p_p , 22,1,1);
      book(_h_all_p_x , 23,1,1);
           
      book(_h_bot_ch_p, 24,1,1);
      book(_h_bot_ch_x, 25,1,1);
      book(_h_bot_pi_p, 26,1,1);
      book(_h_bot_pi_x, 27,1,1);
      book(_h_bot_K_p , 28,1,1);
      book(_h_bot_k_x , 29,1,1);
      book(_h_bot_p_p , 30,1,1);
      book(_h_bot_p_x , 31,1,1);
           
      book(_h_lgt_ch_p, 32,1,1);
      book(_h_lgt_ch_x, 33,1,1);
      book(_h_lgt_pi_p, 34,1,1);
      book(_h_lgt_pi_x, 35,1,1);
      book(_h_lgt_K_p , 36,1,1);
      book(_h_lgt_k_x , 37,1,1);
      book(_h_lgt_p_p , 38,1,1);
      book(_h_lgt_p_x , 39,1,1);

      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<5;++iy) {
	  std::ostringstream title;
	  title << "/TMP/MULT_" << ix << "_" << iy;
	  book(_mult[ix][iy],title.str());
	}
      }
      book(_wLgt,"TMP/wLgt");
      book(_wBot,"TMP/wBot");
      book(_wAll,"TMP/wAll"); 
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<ChargedFinalState>(event, "FS");
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

      // Get event weight for histo filling
      _wAll->fill();
      if(flavour<=3)      _wLgt->fill();
      else if(flavour==5) _wBot->fill();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // loop over the charged particles
      for (const Particle& p : fs.particles()) {
	double modp = p.p3().mod();
        double xp = modp/meanBeamMom;
	int id = abs(p.pid());
	_d_all->fill(modp);
	_mult[0][0]->fill();
	_h_all_ch_p->fill(modp);
	_h_all_ch_x->fill(xp  );
	if(flavour<=3) {
	  _d_lgt->fill(modp);
	  _mult[2][0]->fill();
	  _h_lgt_ch_p->fill(modp);
	  _h_lgt_ch_x->fill(xp  );
	}
	else if(flavour==5) {
	  _d_bot  ->fill(modp);
	  _mult[1][0]->fill();
	  _h_bot_ch_p->fill(modp);
	  _h_bot_ch_x->fill(xp  );
	}
	if(id==211) {
	  _h_all_pi ->fill(modp);
	  _mult[0][1]->fill();
	  _h_all_pi_p->fill(modp);
	  _h_all_pi_x->fill(xp  );
	  if(flavour<=3) {
	    _h_lgt_pi ->fill(modp); 
	    _mult[2][1]->fill();
	    _h_lgt_pi_p->fill(modp);
	    _h_lgt_pi_x->fill(xp  );
	  }
	  else if(flavour==5) {
	    _h_bot_pi ->fill(modp);
	    _mult[1][1]->fill();
	    _h_bot_pi_p->fill(modp);
	    _h_bot_pi_x->fill(xp  );
	  }
	}
	else if(id==321) {
	  _h_all_K ->fill(modp);
	  _h_all_Kp->fill(modp);
	  _mult[0][2]->fill();
	  _mult[0][4]->fill();
	  _h_all_K_p ->fill(modp);
	  _h_all_k_x ->fill(xp  );
	  if(flavour<=3) {
	    _h_lgt_K->fill(modp);
	    _h_lgt_Kp->fill(modp);
	    _mult[2][2]->fill();
	    _mult[2][4]->fill();
	    _h_lgt_K_p ->fill(modp);
	    _h_lgt_k_x ->fill(xp  );
	  }
	  else if(flavour==5) {
	    _h_bot_K ->fill(modp);
	    _h_bot_Kp->fill(modp);
	    _mult[1][2]->fill();
	    _mult[1][4]->fill();
	    _h_bot_K_p ->fill(modp);
	    _h_bot_k_x ->fill(xp  );
	  }
	}
	else if(id==2212) {
	  _h_all_p ->fill(modp);
	  _h_all_Kp->fill(modp);
	  _mult[0][3]->fill();
	  _mult[0][4]->fill();
	  _h_all_p_p ->fill(modp);
	  _h_all_p_x ->fill(xp  );
	  if(flavour<=3) {
	    _h_lgt_p ->fill(modp);
	    _h_lgt_Kp->fill(modp);
	    _mult[2][3]->fill();
	    _mult[2][4]->fill();
	    _h_lgt_p_p ->fill(modp);
	    _h_lgt_p_x ->fill(xp  );
	  }
	  else if(flavour==5) {
	    _h_bot_p ->fill(modp);
	    _h_bot_Kp->fill(modp); 
	    _mult[1][3]->fill();
	    _mult[1][4]->fill();
	    _h_bot_p_p ->fill(modp);
	    _h_bot_p_x ->fill(xp  );
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {


      // // Book histograms
      scale(_h_all_pi,100.);
      scale(_h_all_K ,100.);
      scale(_h_all_p ,100.);
      scale(_h_all_Kp,100.);
      Scatter2DPtr temp;
      book(temp,4,1,1);
      divide(_h_all_pi, _d_all, temp);
      book(temp,5,1,1);
      divide(_h_all_K , _d_all, temp);
      book(temp,6,1,1);
      divide(_h_all_p , _d_all, temp);
      book(temp,7,1,1);
      divide(_h_all_Kp, _d_all, temp);
      
      scale(_h_bot_pi,100.);
      scale(_h_bot_K ,100.);
      scale(_h_bot_p ,100.);
      scale(_h_bot_Kp,100.);
      book(temp, 8,1,1);
      divide(_h_bot_pi, _d_bot, temp);
      book(temp, 9,1,1);
      divide(_h_bot_K , _d_bot, temp);
      book(temp,10,1,1);
      divide(_h_bot_p , _d_bot, temp);
      book(temp,11,1,1);
      divide(_h_bot_Kp, _d_bot, temp);
      
      scale(_h_lgt_pi,100.);
      scale(_h_lgt_K ,100.);
      scale(_h_lgt_p ,100.);
      scale(_h_lgt_Kp,100.);
      book(temp,12,1,1);
      divide(_h_lgt_pi, _d_lgt, temp);
      book(temp,13,1,1);
      divide(_h_lgt_K , _d_lgt, temp);
      book(temp,14,1,1);
      divide(_h_lgt_p , _d_lgt, temp);
      book(temp,15,1,1);
      divide(_h_lgt_Kp, _d_lgt, temp);

      scale(_h_all_ch_p, 1./ *_wAll);
      scale(_h_all_ch_x, 1./ *_wAll);
      scale(_h_all_pi_p, 1./ *_wAll);
      scale(_h_all_pi_x, 1./ *_wAll);
      scale(_h_all_K_p , 1./ *_wAll);
      scale(_h_all_k_x , 1./ *_wAll);
      scale(_h_all_p_p , 1./ *_wAll);
      scale(_h_all_p_x , 1./ *_wAll);

      scale(_h_bot_ch_p, 1./ *_wBot);
      scale(_h_bot_ch_x, 1./ *_wBot);
      scale(_h_bot_pi_p, 1./ *_wBot);
      scale(_h_bot_pi_x, 1./ *_wBot);
      scale(_h_bot_K_p , 1./ *_wBot);
      scale(_h_bot_k_x , 1./ *_wBot);
      scale(_h_bot_p_p , 1./ *_wBot);
      scale(_h_bot_p_x , 1./ *_wBot);

      scale(_h_lgt_ch_p, 1./ *_wLgt);
      scale(_h_lgt_ch_x, 1./ *_wLgt);
      scale(_h_lgt_pi_p, 1./ *_wLgt);
      scale(_h_lgt_pi_x, 1./ *_wLgt);
      scale(_h_lgt_K_p , 1./ *_wLgt);
      scale(_h_lgt_k_x , 1./ *_wLgt);
      scale(_h_lgt_p_p , 1./ *_wLgt);
      scale(_h_lgt_p_x , 1./ *_wLgt);

      // multiplicities
      vector<CounterPtr> scales = {_wAll,_wBot,_wLgt};
      for(unsigned int ix=0;ix<3;++ix) {
	if(scales[ix]->effNumEntries()<=0.) continue;
	for(unsigned int iy=0;iy<5;++iy) {
	  Scatter2DPtr scatter;
	  book(scatter, ix+1, 1, iy+1, true);
	  scale(_mult[ix][iy],1./ *scales[ix]);
	  scatter->point(0).setY(_mult[ix][iy]->val(),_mult[ix][iy]->err());
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_all_pi , _h_all_K  , _h_all_p  , _h_all_Kp , _d_all;
    Histo1DPtr _h_bot_pi , _h_bot_K  , _h_bot_p  , _h_bot_Kp , _d_bot;
    Histo1DPtr _h_lgt_pi , _h_lgt_K  , _h_lgt_p  , _h_lgt_Kp , _d_lgt;
    Histo1DPtr _h_all_ch_p, _h_all_ch_x , _h_all_pi_p , _h_all_pi_x ;
    Histo1DPtr _h_all_K_p , _h_all_k_x  , _h_all_p_p  , _h_all_p_x  ;
    Histo1DPtr _h_bot_ch_p , _h_bot_ch_x , _h_bot_pi_p , _h_bot_pi_x;
    Histo1DPtr _h_bot_K_p  , _h_bot_k_x  , _h_bot_p_p  , _h_bot_p_x ;
    Histo1DPtr _h_lgt_ch_p , _h_lgt_ch_x , _h_lgt_pi_p , _h_lgt_pi_x;
    Histo1DPtr _h_lgt_K_p  , _h_lgt_k_x  , _h_lgt_p_p  , _h_lgt_p_x ;
    CounterPtr _mult[3][5];

    CounterPtr _wLgt, _wBot, _wAll;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(DELPHI_1998_I473409);


}
