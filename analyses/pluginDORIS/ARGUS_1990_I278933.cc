// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0 and eta at Upsilon 1,2 and continuum
  class ARGUS_1990_I278933 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1990_I278933);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_cont_pi1 ,3, 1, 1);
      book(_h_cont_pi2 ,3, 1, 2);
      book(_h_ups1_pi  ,4, 1, 1);
      book(_h_ups2_pi  ,4, 1, 2);
      book(_h_cont_eta1,5, 1, 1);
      book(_h_cont_eta2,5, 1, 2);
      book(_h_ups1_eta ,6, 1, 1);
      book(_h_ups2_eta ,6, 1, 2);
      book(_n_Eta[0], "/TMP/EtaCont");
      book(_n_Eta[1], "/TMP/EtaUps1");
      book(_n_Eta[2], "/TMP/EtaUps2");
      book(_n_Pi[0] , "/TMP/PiCont");
      book(_n_Pi[1] , "/TMP/PiUps1");
      book(_n_Pi[2] , "/TMP/PiUps2");
      book(_weightSum_cont,"/TMP/weightSum_cont");
      book(_weightSum_Ups1,"/TMP/weightSum_Ups1");
      book(_weightSum_Ups2,"/TMP/weightSum_Ups2");
		      
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
	if(id == 111 or id == 221) {
	  unstable.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 or Cuts::pid==100553);
      // Continuum
      if (upsilons.empty()) { 
        MSG_DEBUG("No Upsilons found => continuum event");
        _weightSum_cont->fill();
        for (const Particle& p : ufs.particles(Cuts::pid==111 or Cuts::pid==221)) {
          const int id = p.pid();
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
	  if(id==111) {
	    _n_Pi[0]->fill();
	    _h_cont_pi1->fill(xp,1./beta);
	    _h_cont_pi2->fill(xp,1./beta);
	  }
	  else {
	    _n_Eta[0]->fill();
	    _h_cont_eta1->fill(xp,1./beta);
	    _h_cont_eta2->fill(xp,1./beta);
	  }
	}
      }
      // Upsilon(s) found
      else { 
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
          const int parentId = ups.pid();
	  if(parentId==553) {
	    _weightSum_Ups1->fill();
	  }
	  else {
	    _weightSum_Ups2->fill();
	  }
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups, unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
	  // loop over decay products
          for(const Particle& p : unstable) {
            const int id = p.pid();
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
	    if(id==111) {
	      if(parentId==553) {
		_n_Pi[1]->fill();
		_h_ups1_pi->fill(xp,1./beta);
	      }
	      else {
		_n_Pi[2]->fill();
		_h_ups2_pi->fill(xp,1./beta);
	      }
	    }
	    else if(id==221) {
	      if(parentId==553) {
		_n_Eta[1]->fill();
		_h_ups1_eta->fill(xp,1./beta);
	      }
	      else {
		_n_Eta[2]->fill();
		_h_ups2_eta->fill(xp,1./beta);
	      }
	    } 
	  }
	}
      }	
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Scale histos
      if (_weightSum_cont->val() > 0.) {
	scale(_h_cont_pi1, 1./ *_weightSum_cont);
	scale(_h_cont_pi2, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
	scale(_h_cont_eta1, 1./ *_weightSum_cont);
	scale(_h_cont_eta2, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      }
      if (_weightSum_Ups1->val() > 0.) {
	scale(_h_ups1_pi, 1./ *_weightSum_Ups1);
	scale(_h_ups1_eta, 1./ *_weightSum_Ups1);
      }
      if (_weightSum_Ups2->val() > 0.) {
	scale(_h_ups2_pi, 1./ *_weightSum_Ups2);
	scale(_h_ups2_eta, 1./ *_weightSum_Ups2);
      }
      // Counters
      vector<CounterPtr> scales = {_weightSum_cont,_weightSum_Ups1,_weightSum_Ups2};
      for(unsigned int ix=0;ix<3;++ix) {
	Scatter2DPtr scatter;
	book(scatter, 1, 1, ix+1, true);
	scale(_n_Pi[ix],1./ *scales[ix]);
        scatter->point(0).setY(_n_Pi[ix]->val(),
			       _n_Pi[ix]->err());
      }
      for(unsigned int ix=0;ix<3;++ix) {
	Scatter2DPtr scatter;
	book(scatter, 2, 1, ix+1, true);
	scale(_n_Eta[ix],1./ *scales[ix]);
        scatter->point(0).setY(_n_Eta[ix]->val(),
			       _n_Eta[ix]->err());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_cont_pi1 , _h_cont_pi2 , _h_ups1_pi , _h_ups2_pi ;
    Histo1DPtr _h_cont_eta1, _h_cont_eta2, _h_ups1_eta, _h_ups2_eta;
    CounterPtr _n_Eta[3],_n_Pi[3];
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1990_I278933);


}
