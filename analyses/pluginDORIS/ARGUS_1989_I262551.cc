// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief phi spectrum in continuum, and Upsilon 1s and 2s decays
  class ARGUS_1989_I262551 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1989_I262551);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_cont, 1, 1, 1);
      book(_h_ups1, 2, 1, 1);
      book(_h_ups2, 2, 1, 2);
      book(_n_Phi[0], "/TMP/NUps1");
      book(_n_Phi[1], "/TMP/NUps2");
      book(_weightSum_cont, "TMP/weightSum_cont");
      book(_weightSum_Ups1, "TMP/weightSum_Ups1");
      book(_weightSum_Ups2, "TMP/weightSum_Ups2");
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& phis) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
	if(id == 333) {
	  phis.push_back(p);
	}
	if(!p.children().empty())
	  findDecayProducts(p, phis);
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
        for (const Particle& p : ufs.particles(Cuts::pid==333)) {
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
	  _h_cont->fill(xp,1./beta);
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
          Particles phis;
          // Find the decay products we want
          findDecayProducts(ups, phis);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
	  // loop over decay products
          for(const Particle& p : phis) {
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
	    if(parentId==553) {
	      _n_Phi[0]->fill();
	      _h_ups1->fill(xp,1./beta);
	    }
	    else {
	      _n_Phi[1]->fill();
	      _h_ups2->fill(xp,1./beta);
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_weightSum_cont->effNumEntries() > 0.)
	scale(_h_cont, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      if (_weightSum_Ups1->effNumEntries() > 0.) {
	scale(_h_ups1, 1./ *_weightSum_Ups1);
      }
      if (_weightSum_Ups2->effNumEntries() > 0.) {
	scale(_h_ups2, 1./ *_weightSum_Ups2);
      }
      // Counters
      vector<CounterPtr> scales = {_weightSum_Ups1,_weightSum_Ups2};
      for(unsigned int ix=0;ix<2;++ix) {
	Scatter2DPtr scatter;
	book(scatter, 3+ix, 1, 1, true);
	scale(_n_Phi[ix], 1./ *scales[ix]);
        scatter->point(0).setY(_n_Phi[ix]->val(),
			       _n_Phi[ix]->err());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_cont, _h_ups1, _h_ups2;
    CounterPtr _n_Phi[2];
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ARGUS_1989_I262551);


}
