// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Thrust like variable at Upsilon(1s,2S)
  class LENA_1981_I164397 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LENA_1981_I164397);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      declare(ChargedFinalState(), "FS");
      
      book(_weightSum_cont, "TMP/weightSum_cont");
      book(_weightSum_Ups1, "TMP/weightSum_Ups1");
      book(_weightSum_Ups2, "TMP/weightSum_Ups2");
      book(_charge_cont, "TMP/charge_cont");
      book(_charge_Ups1, "TMP/charge_Ups1");
      book(_charge_Ups2, "TMP/charge_Ups2");

      if(fuzzyEquals(sqrtS(),9.5149,1e-2)) {
	book(_hist_T_cont ,4, 1, 1);
      }
      else if(fuzzyEquals(sqrtS(),9.9903,1e-2)) {
	book(_hist_T_cont ,4, 1, 2);
      }
      book(_hist_T_Ups1 ,4, 1, 3);
      book(_hist_T_Ups2 ,4, 1, 4);
    }

    /// Recursively walk the decay tree to find the charged decay products of @a p
    void findDecayProducts(Particle mother, Particles& charged) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
	if(!p.children().empty())
	  findDecayProducts(p, charged);
	else if(PID::isCharged(id))
	  charged.push_back(p);
      }
    }

    // defn of thrust in paper used just the direction
    double thrustPrime(const LorentzTransform & boost, const Particles & particles) {
      vector<Vector3> vecs;
      for(const Particle & p : particles) {
	vecs.push_back(boost.transform(p.momentum()).p3().unit());
      }
      Thrust thrust;
      thrust.calc(vecs);
      return thrust.thrust();
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 or Cuts::pid==100553);
      if (upsilons.empty()) { 
        MSG_DEBUG("No Upsilons found => continuum event");
        _weightSum_cont->fill();
	Particles cfs = apply<ChargedFinalState>(event, "FS").particles();
	_charge_cont->fill(cfs.size());
	if(_hist_T_cont) {
          LorentzTransform boost;
	  _hist_T_cont->fill(thrustPrime(boost,cfs));
	}
      }
      // Upsilon(s) found
      else {
        for (const Particle& ups : upsilons) {
          const int parentId = ups.pid();
          Particles charged;
	  // boost to rest frame (if required)
          LorentzTransform boost;
          if (ups.p3().mod() > 1*MeV)
            boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          // Find the decay products we want
          findDecayProducts(ups, charged);
	  if(parentId==553) {
	    _weightSum_Ups1->fill();
	    _charge_Ups1->fill(charged.size());
	    _hist_T_Ups1->fill(thrustPrime(boost,charged));
	  }
	  else {
	    _weightSum_Ups2->fill();
	    _charge_Ups2->fill(charged.size());
	    _hist_T_Ups2->fill(thrustPrime(boost,charged));
	  }
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // charged particle multiplicity
      if(_weightSum_cont->val()>0. ) {
	scale(_charge_cont,1./ *_weightSum_cont );
	if(_hist_T_cont) scale(_hist_T_cont,1./ *_weightSum_cont );
      }
      if(_weightSum_Ups1->val()>0. ) {
	scale(_charge_Ups1,1./ *_weightSum_Ups1 );
	 scale(_hist_T_Ups1,1./ *_weightSum_Ups1 );
      }
      if(_weightSum_Ups2->val()>0. ) {
	scale(_charge_Ups2,1./ *_weightSum_Ups2 );
	scale(_hist_T_Ups2,1./ *_weightSum_Ups2 );
      }
      Scatter2D tempScat(refData(3, 1, 1));
      Scatter2DPtr _mult;
      book(_mult, 3, 1, 1);
      for (size_t b = 0; b < tempScat.numPoints(); b++) {
        const double x  = tempScat.point(b).x();
        pair<double,double> ex = tempScat.point(b).xErrs();
        pair<double,double> ex2 = ex;
        if(ex2.first ==0.) ex2. first=0.02;
        if(ex2.second==0.) ex2.second=0.02;
	// Upsilon 1S
	if(b==3) {
	  if (_weightSum_Ups1->val()>0.) {
	    _mult->addPoint(x, _charge_Ups1->val(), ex, make_pair(_charge_Ups1->err(),_charge_Ups1->err()));
	  }
	  else {
	    _mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
	// Upsilon 2S
	else if(b==6) {
	  if (_weightSum_Ups2->val()>0.) {
	    _mult->addPoint(x, _charge_Ups2->val(), ex, make_pair(_charge_Ups2->err(),_charge_Ups2->err()));
	  }
	  else {
	    _mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
	else {
	  if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second) && _weightSum_cont->val()>0.) {
	    _mult->addPoint(x, _charge_cont->val(), ex, make_pair(_charge_cont->err(),_charge_cont->err()));
	  }
	  else {
	    _mult->addPoint(x, 0., ex, make_pair(0.,.0));
	  }
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    CounterPtr _weightSum_cont, _weightSum_Ups1, _weightSum_Ups2;
    CounterPtr _charge_cont, _charge_Ups1, _charge_Ups2;
    Histo1DPtr _hist_T_cont,_hist_T_Ups1,_hist_T_Ups2;
    //@}


  };


  DECLARE_RIVET_PLUGIN(LENA_1981_I164397);

}
