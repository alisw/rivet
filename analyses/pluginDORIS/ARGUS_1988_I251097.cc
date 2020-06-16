// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief hyperon production
  class ARGUS_1988_I251097 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1988_I251097);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<7;++iy) {
	  std::ostringstream title;
	  title << "/TMP/MULT_" << ix << "_" << iy;
	  book(_mult[ix][iy],title.str());
	}
      }
      book(_hist_ups1_lambda , 3,1,1);
      book(_hist_ups2_lambda , 4,1,1);
      book(_hist_cont_lambda1, 5,1,1);
      book(_hist_cont_lambda2, 6,1,1);
           
      book(_hist_ups1_xi     , 7,1,1);
      book(_hist_ups2_xi     , 8,1,1);
      book(_hist_cont_xi     , 9,1,1);
      book(_weightSum_cont,"TMP/sum_cont");
      book(_weightSum_Ups1,"TMP/sum_ups1");
      book(_weightSum_Ups2,"TMP/sum_ups2");
    }

   /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
	if (id == 3122 || id == 3312 || id == 3212 || id == 3114 ||
            id == 3224 || id == 3324 || id == 3334) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First in unstable final state
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 || Cuts::pid==100553);
      // continuum
      if (upsilons.empty()) { 
        _weightSum_cont->fill();
        for (const Particle& p : ufs.particles()) {
          int id = p.abspid();
	  double modp = p.p3().mod();
          double xp = 2.*modp/sqrtS();
          double xE = 2.*p.E()/sqrtS();
          double beta = modp/p.E();
          if (id == 3122) {
	    _hist_cont_lambda1->fill(xp);
	    _hist_cont_lambda2->fill(xE, 1./beta);
	    _mult[1][0]->fill();
          }
          else if (id == 3312) {
            _hist_cont_xi->fill(xE, 1./beta);
	    _mult[1][1]->fill();
          }
          else if (id == 3212) {
	    _mult[1][2]->fill();
          }
          else if (id == 3114) {
	    _mult[1][3]->fill();
          }
          else if (id == 3224) {
	    _mult[1][4]->fill();
          }
          else if (id == 3324) {
	    _mult[1][5]->fill();
          }
          else if (id == 3334) {
	    _mult[1][6]->fill();
          }
        }
      }
      // Upslion decays
      else {
	for (const Particle& ups : upsilons) {
	  const int parentId = ups.pid();
	  if( parentId == 553)
	    _weightSum_Ups1->fill();
	  else
	    _weightSum_Ups2->fill();
	  Particles unstable;
	  // Find the decay products we want
	  findDecayProducts(ups,unstable);
	  LorentzTransform cms_boost;
	  if (ups.p3().mod() > 0.001)
	    cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	  double mass = ups.mass();
	  for(const Particle & p : unstable) {
	    int id = p.abspid();
	    FourMomentum p2 = cms_boost.transform(p.momentum());
	    double modp = p2.p3().mod();
	    double xp = 2.*modp/mass;
	    if (id == 3122) {
	      if(parentId==553) {
		_hist_ups1_lambda->fill(xp);
		_mult[0][0]->fill();
	      }
	      else {
		_hist_ups2_lambda->fill(xp);
	      }
	    }
	    else if (id == 3312) {
	      if(parentId==553) {
		_hist_ups1_xi->fill(xp);
		_mult[0][1]->fill();
	      }
	      else {
		_hist_ups2_xi->fill(xp);
	      }
	    }
	    else if(parentId==553) {
	      if (id == 3212) {
		_mult[0][2]->fill();
	      }
	      else if (id == 3114) {
		_mult[0][3]->fill();
	      }
	      else if (id == 3224) {
		_mult[0][4]->fill();
	      }
	      else if (id == 3324) {
		_mult[0][5]->fill();
	      }
	      else if (id == 3334) {
		_mult[0][6]->fill();
	      }
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // multiplicities
      vector<CounterPtr> scales = {_weightSum_Ups1,_weightSum_cont};
      for(unsigned int ix=0;ix<2;++ix) {
	if(scales[ix]->effNumEntries()<=0.) continue;
	for(unsigned int iy=0;iy<7;++iy) {
	  Scatter2DPtr scatter;
	  book(scatter, ix+1, 1, iy+1, true);
	  scale(_mult[ix][iy],1./ *scales[ix]);
	  scatter->point(0).setY(_mult[ix][iy]->val(),_mult[ix][iy]->err());
	}
      }
      if(_weightSum_Ups1->val()>0.) {
	scale(_hist_ups1_lambda,1./ *_weightSum_Ups1);
	scale(_hist_ups1_xi    ,1./ *_weightSum_Ups1);
      }
      if(_weightSum_Ups2->val()>0.) {
	scale(_hist_ups2_lambda,1./ *_weightSum_Ups2);
	scale(_hist_ups2_xi    ,1./ *_weightSum_Ups2);
      }
      if(_weightSum_cont->val()) {
	scale(_hist_cont_lambda1, 1./ *_weightSum_cont);
	scale(_hist_cont_lambda2, 1./ *_weightSum_cont);
	scale(_hist_cont_xi     , 1./ *_weightSum_cont);
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist_ups1_lambda, _hist_ups2_lambda, _hist_cont_lambda1, _hist_cont_lambda2;
    Histo1DPtr _hist_ups1_xi, _hist_ups2_xi, _hist_cont_xi;
    CounterPtr _mult[2][7];
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1988_I251097);


}
