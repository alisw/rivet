// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0 and eta at Upsilon 1 and continuum
  class CRYSTAL_BALL_1991_I297905 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1991_I297905);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_cont_pi[0] ,3, 1, 1);
      book(_h_cont_pi[1] ,3, 1, 2);
      book(_h_ups1_pi    ,4, 1, 1);
      book(_h_cont_eta[0],5, 1, 1);
      book(_h_cont_eta[1],5, 1, 2);
      book(_h_ups1_eta ,6, 1, 1);
      // counters
      book(_n_Eta[0], "/TMP/EtaCont");
      book(_n_Eta[1], "/TMP/EtaUps1");
      book(_n_Pi[0] , "/TMP/PiCont");
      book(_n_Pi[1] , "/TMP/PiUps1");
      book(_weightSum_cont,"/TMP/weightSum_cont");
      book(_weightSum_Ups1,"/TMP/weightSum_Ups1");
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
      Particles upsilons = ufs.particles(Cuts::pid==553);
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
	    for(unsigned int ix=0;ix<2;++ix)
	      _h_cont_pi[ix]->fill(xp,1./beta);
	  }
	  else {
	    _n_Eta[0]->fill();
	    for(unsigned int ix=0;ix<2;++ix)
	      _h_cont_eta[ix]->fill(xp,1./beta);
	  }
	}
      }
      // Upsilon(s) found
      else { 
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
	  _weightSum_Ups1->fill();
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
	      _n_Pi[1]->fill();
	      _h_ups1_pi->fill(xp,1./beta);
	    }
	    else if(id==221) {
	      _n_Eta[1]->fill();
	      _h_ups1_eta->fill(xp,1./beta);
	    } 
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Scale histos
      if (_weightSum_cont->val() > 0.) {
	scale(_h_cont_pi[0] , 1./ *_weightSum_cont);
	scale(_h_cont_pi[1] , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
	scale(_h_cont_eta[0], 1./ *_weightSum_cont);
	scale(_h_cont_eta[1], sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
	Scatter2DPtr scatter;
       	book(scatter, 1, 1, 1, true);
      	scale(_n_Pi[0],1./ *_weightSum_cont);
        scatter->point(0).setY(_n_Pi[0]->val(),
      			       _n_Pi[0]->err());
       	book(scatter, 1, 1, 2, true);
      	scale(_n_Eta[0],1./ *_weightSum_cont);
        scatter->point(0).setY(_n_Eta[0]->val(),
      			       _n_Eta[0]->err());
      }
      if (_weightSum_Ups1->val() > 0.) {
	scale(_h_ups1_pi, 1./ *_weightSum_Ups1);
	scale(_h_ups1_eta, 1./ *_weightSum_Ups1);
	Scatter2DPtr scatter;
       	book(scatter, 2, 1, 1, true);
      	scale(_n_Pi[1],1./ *_weightSum_Ups1);
        scatter->point(0).setY(_n_Pi[1]->val(),
      			       _n_Pi[1]->err());
       	book(scatter, 2, 1, 2, true);
      	scale(_n_Eta[1],1./ *_weightSum_Ups1);
        scatter->point(0).setY(_n_Eta[1]->val(),
      			       _n_Eta[1]->err());
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cont_pi[2] , _h_ups1_pi;
    Histo1DPtr _h_cont_eta[2], _h_ups1_eta;
    CounterPtr _n_Eta[2],_n_Pi[2];
    CounterPtr _weightSum_cont,_weightSum_Ups1;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1991_I297905);

}
