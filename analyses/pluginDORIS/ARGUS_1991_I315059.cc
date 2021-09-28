// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief ARGUS D0, D+, D*+ production
  class ARGUS_1991_I315059 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1991_I315059);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      declare(Beam(), "Beams");

      // Book histograms
      book(_n_D0   , 1,1,1);
      book(_n_Dp   , 1,1,2);
      book(_n_DStar, 1,1,3);

      book(_h_x_D0   , 2,1,1);
      book(_h_x_Dp   , 3,1,1);
      book(_h_x_DStar, 4,1,1);

      book(_h_p_D0   , 5,1,1);
      book(_h_p_Dp   , 6,1,1);
      book(_h_p_DStar, 7,1,1);

      book(_c_cont, "/TMP/c_cont");
      book(_c_ups , "/TMP/c_ups");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // Find the Upsilon(4S) among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==300553);
      // Continuum
      if (upsilons.empty()) {
	_c_cont->fill();
        MSG_DEBUG("No Upsilons found => continuum event");
	for(const Particle &p : ufs.particles(Cuts::abspid==411 or Cuts::abspid==421 or Cuts::abspid==413)) {
	  double xp = p.p3().mod()/sqrt(sqr(meanBeamMom)+sqr(p.mass()));
	  if(p.abspid()==421) {
	    _h_x_D0->fill(xp);
	    _n_D0->fill(10.6);
	  }
	  else if(p.abspid()==411) {
	    _h_x_Dp->fill(xp);
	    _n_Dp->fill(10.6);
	  }
	  else if(p.abspid()==413) {
	    _h_x_DStar->fill(xp);
	    _n_DStar->fill(10.6);
	  }
	}
      }
      // upsilon decay
      else {
        for (const Particle& ups : upsilons) {
	  _c_ups->fill();
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups, unstable);
	  // boost to rest frame (if required)
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          for(const Particle& p : unstable) {
            const FourMomentum p2 = cms_boost.transform(p.momentum());
	    double modp = p2.p3().mod();
	    if(p.abspid()==421) {
	      _h_p_D0->fill(modp);
	    }
	    else if(p.abspid()==411) {
	      _h_p_Dp->fill(modp);
	    }
	    else if(p.abspid()==413) {
	      _h_p_DStar->fill(modp);
	    }
	  }
	}
      }
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle & p: mother.children()) {
        const int id = p.abspid();
	if (id == 411 || id == 421) {
	  unstable.push_back(p);
	}
	else if (id == 413 ) {
	  unstable.push_back(p);
	  findDecayProducts(p, unstable);
	}
	else if(!p.children().empty())
	  findDecayProducts(p, unstable);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // brs for the decays used (PDG 2018)
      double brD0 = 0.0389;
      double brDp = 0.0898;
      double brDStar = 0.677;
      if(_c_cont->numEntries()!=0) {
	scale(_n_D0      , 1./sumOfWeights()*crossSection()/nanobarn);
	scale(_n_Dp      , 1./sumOfWeights()*crossSection()/nanobarn);
	scale(_n_DStar  ,  1./sumOfWeights()*crossSection()/nanobarn);
	scale(_h_x_D0   , brD0/sumOfWeights()*crossSection()/nanobarn*sqr(sqrtS()));
	scale(_h_x_Dp   , brDp/sumOfWeights()*crossSection()/nanobarn*sqr(sqrtS()));
	scale(_h_x_DStar, brD0*brDStar/sumOfWeights()*crossSection()/nanobarn*sqr(sqrtS()));
      }
      if(_c_ups->numEntries()!=0) {
	scale(_h_p_D0   , 1e3*0.5*brD0/ *_c_ups);
	scale(_h_p_Dp   , 1e3*0.5*brDp/ *_c_ups);
	scale(_h_p_DStar, 1e3*0.5*brD0*brDStar/ *_c_ups);
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _n_D0,_n_Dp,_n_DStar;
    Histo1DPtr _h_x_D0,_h_x_Dp,_h_x_DStar;
    Histo1DPtr _h_p_D0,_h_p_Dp,_h_p_DStar;
    CounterPtr _c_cont, _c_ups;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1991_I315059);


}
