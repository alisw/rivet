// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Production of the $\eta'(958)$ and $f_0(980)$ in $e^+e^-$ annihilation in the Upsilon region
  ///
  /// @author Peter Richardson
  class ARGUS_1993_S2669951 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1993_S2669951);


    void init() {
      declare(UnstableParticles(), "UFS");

      book(_weightSum_cont, "TMP/weightSum_cont");
      book(_weightSum_Ups1, "TMP/weightSum_Ups1");
      book(_weightSum_Ups2, "TMP/weightSum_Ups2");

      for ( auto i : {0,1,2} ) {
        if ( i < 2 )
          book(_count_etaPrime_highZ[i], "TMP/count_etaPrime_highz_" + to_str(i));
        book(_count_etaPrime_allZ[i], "TMP/count_etaPrime_allz_" + to_str(i));
        book(_count_f0[i], "TMP/count_f0_" + to_str(i));
      }

      book(_hist_cont_f0 ,2, 1, 1);
      book(_hist_Ups1_f0 ,3, 1, 1);
      book(_hist_Ups2_f0 ,4, 1, 1);
    }


    void analyze(const Event& e) {
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(e, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553 or Cuts::pid==100553);
      // Continuum
      if (upsilons.empty()) {
        MSG_DEBUG("No Upsilons found => continuum event");
        _weightSum_cont->fill();
        for (const Particle& p : ufs.particles()) {
          const int id = p.pid();
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
          if (id == 9010221) {
            _hist_cont_f0->fill(xp, 1./beta);
	    _count_f0[2]->fill();
          } else if (id == 331) {
            if (xp > 0.35) _count_etaPrime_highZ[1]->fill();
	    _count_etaPrime_allZ[2]->fill();
          }
	}
      }
      // Upsilon(s) found
      else {
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
	  // boost to rest frame (if required)
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
            if (id == 9010221) {
	      if(parentId == 553 ) {
		_hist_Ups1_f0->fill(xp, 1./beta);
		_count_f0[0]->fill();
	      }
	      else {
		_hist_Ups2_f0->fill(xp, 1./beta);
		_count_f0[1]->fill();
	      }
	    }
	    else if ( id == 331 ) {
	      if (parentId == 553) {
		if (xp > 0.35) _count_etaPrime_highZ[0]->fill();
		_count_etaPrime_allZ[0]->fill();
	      }
	      else {
		_count_etaPrime_allZ[1]->fill();
	      }
	    }
	  }
	}
      }
    }

    void finalize() {
      // High-Z eta' multiplicity
      Scatter2DPtr s111;
      book(s111, 1, 1, 1, true);
      if (_weightSum_Ups1->val() > 0) // Point at 9.460
        s111->point(0).setY(_count_etaPrime_highZ[0]->val() / _weightSum_Ups1->val(), 0);
      if (_weightSum_cont->val() > 0) // Point at 9.905
        s111->point(1).setY(_count_etaPrime_highZ[1]->val() / _weightSum_cont->val(), 0);

      // All-Z eta' multiplicity
      Scatter2DPtr s112;
      book(s112, 1, 1, 2, true);
      if (_weightSum_Ups1->val() > 0) // Point at 9.460
        s112->point(0).setY(_count_etaPrime_allZ[0]->val() / _weightSum_Ups1->val(), 0);
      if (_weightSum_cont->val() > 0) // Point at 9.905
        s112->point(1).setY(_count_etaPrime_allZ[2]->val() / _weightSum_cont->val(), 0);
      if (_weightSum_Ups2->val() > 0) // Point at 10.02
        s112->point(2).setY(_count_etaPrime_allZ[1]->val() / _weightSum_Ups2->val(), 0);


      // f0 multiplicity
      Scatter2DPtr s511;
      book(s511, 5, 1, 1, true);
      if (_weightSum_Ups1->val() > 0) // Point at 9.46
        s511->point(0).setY(_count_f0[0]->val() / _weightSum_Ups1->val(), 0);
      if (_weightSum_Ups2->val() > 0) // Point at 10.02
        s511->point(1).setY(_count_f0[1]->val() / _weightSum_Ups2->val(), 0);
      if (_weightSum_cont->val() > 0) // Point at 10.45
        s511->point(2).setY(_count_f0[2]->val() / _weightSum_cont->val(), 0);

      // Scale histos
      if (_weightSum_cont->val() > 0.) scale(_hist_cont_f0, 1./ *_weightSum_cont);
      if (_weightSum_Ups1->val() > 0.) scale(_hist_Ups1_f0, 1./ *_weightSum_Ups1);
      if (_weightSum_Ups2->val() > 0.) scale(_hist_Ups2_f0, 1./ *_weightSum_Ups2);
    }


  private:

    /// @name Counters
    /// @{
    array<CounterPtr,3> _count_etaPrime_highZ, _count_etaPrime_allZ, _count_f0;
    CounterPtr _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    /// @}


    /// Histos
    Histo1DPtr _hist_cont_f0, _hist_Ups1_f0, _hist_Ups2_f0;


    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
	if (id == 331 || id == 9010221) {
	  unstable.push_back(p);
	}
	else if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

  };



  RIVET_DECLARE_ALIASED_PLUGIN(ARGUS_1993_S2669951, ARGUS_1993_I342061);

}
