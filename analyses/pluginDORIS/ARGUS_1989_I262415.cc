// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda(1520) production
  class ARGUS_1989_I262415 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1989_I262415);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_ups1_obs, 1, 1, 1);
      book(_h_cont_obs, 1, 1, 2);
      book(_h_ups1_all, 1, 2, 1);
      book(_h_cont_all, 1, 2, 2);
      book(_h_ups1    , 3, 1, 1);
      book(_h_cont    , 4, 1, 1);
      book(_weightSum_cont,"TMP/weightSum_cont");
      book(_weightSum_Ups1,"TMP/weightSum_Ups1");

    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
	if(id == 3124) {
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
        for (const Particle& p : ufs.particles(Cuts::abspid==3124)) {
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
	  _h_cont->fill(xp,1./beta);
	  _h_cont_obs->fill(xp);
	  _h_cont_all->fill(xp);
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
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
	    _h_ups1->fill(xp,1./beta);
	    _h_ups1_obs->fill(xp);
	    _h_ups1_all->fill(xp);
	  }
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      if (_weightSum_cont->val() > 0.) {
	scale(_h_cont    , 1. / *_weightSum_cont);
	scale(_h_cont_obs, 0.3/ *_weightSum_cont);
	scale(_h_cont_all, 1. / *_weightSum_cont);
      }
      if (_weightSum_Ups1->val() > 0.) {
	scale(_h_ups1    , 1. / *_weightSum_Ups1);
	scale(_h_ups1_obs, 0.3/ *_weightSum_Ups1);
	scale(_h_ups1_all, 1. / *_weightSum_Ups1);
      }

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ups1, _h_cont;
    Histo1DPtr _h_ups1_obs, _h_cont_obs;
    Histo1DPtr _h_ups1_all, _h_cont_all;
    CounterPtr _weightSum_cont,_weightSum_Ups1;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1989_I262415);


}
