// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta' production in Upsilon(1S) decays
  class CLEOIII_2006_I728679 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOIII_2006_I728679);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");

      book(_hist    , 1, 1, 1);
      book(_mult,"TMP/mult");
      book(_weightSum,"TMP/weightSum");
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = abs(p.pid());
	if ( id == 331 ) {
	  unstable.push_back(p);
	}
	else if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the Upsilon(1S) mesons
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      if(upsilons.empty()) vetoEvent;
      // loop over them
      for (const Particle& ups : upsilons) {
	_weightSum->fill();
	Particles unstable;
	// Find the decay products we want
	findDecayProducts(ups,unstable);
	LorentzTransform cms_boost;
	if (ups.p3().mod() > 0.001)
	  cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
	double mass = ups.mass();
	_mult->fill(unstable.size());
	// fill the spectrum
	for( const Particle & p : unstable) {
	  FourMomentum p2 = cms_boost.transform(p.momentum());
	  double xp = 2.*p2.E()/mass;
	  _hist->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // spectrum
      if (_weightSum->val() > 0.)
        scale(_hist, 1. / *_weightSum);
      // BR
      Scatter2DPtr scatter;
      book(scatter,2, 1, 1, true);
      scale(_mult,1./ *_weightSum);
      scatter->point(0).setY(_mult->val(),_mult->err());
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _hist;
    CounterPtr _mult;
    CounterPtr _weightSum;
    //@}


  };


  DECLARE_RIVET_PLUGIN(CLEOIII_2006_I728679);

}
