// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class HRS_1990_I280958 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(HRS_1990_I280958);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Thrust(cfs), "Thrust");

      // Book histograms
      book(_h_X        , 3, 1, 1);
      book(_h_rap_all  , 4, 1, 1);
      book(_h_rap_light, 6, 1, 1);
      book(_h_rap_charm, 5, 1, 1);
      book(_wLight,"TMP/wLight");
      book(_wCharm,"TMP/wCharm");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      int nch = cfs.particles().size();
      if(nch<5) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // get the thrust axes
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const Vector3 & axis = thrust.thrustAxis();
      // unstable particles
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particle pTag;
      // get the tags
      Particles Dstar = ufs.particles(Cuts::abspid==413);
      bool charmTagged = !Dstar.empty();
      if(charmTagged) {
	pTag = Dstar[0];
	for(const Particle & p : Dstar) {
	  if(p.E()>pTag.E()) pTag=p;
	}
	_wCharm->fill();
      }
      bool lightTagged = false;
      if(!charmTagged) {
	for(const Particle & p : cfs.particles()) {
	  if(p.p3().mod()>9.43*GeV) {
	    pTag=p;
	    lightTagged=true;
	    _wLight->fill();
	    break;
	  }
	}
      }
      // sign of hemispheres if tagged
      double sign=1.;
      if(charmTagged || lightTagged) {
	if(dot(axis,pTag.p3())<0.) sign=-1.;
      }
      // now loop over the kaons
      for(const Particle & p : ufs.particles(Cuts::pid==130 || Cuts::pid==310)) {
         double xE = p.E()/meanBeamMom;
	 const double energy = p.E();
	 const double momT = dot(axis, p.p3());
	 _h_X->fill(xE);
	 double rap = 0.5 * std::log((energy + momT) / (energy - momT));
	 _h_rap_all->fill(fabs(rap));
	 rap *=sign;
	 if(charmTagged)
	   _h_rap_charm->fill(rap);
	 else if(lightTagged)
	   _h_rap_light->fill(rap);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_X        , crossSection()*sqr(sqrtS())/nanobarn/sumOfWeights());
      scale(_h_rap_all  , 1./sumOfWeights());
      scale(_h_rap_light, 1./ *_wLight);
      scale(_h_rap_charm, 1./ *_wCharm);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_X, _h_rap_all,_h_rap_light,_h_rap_charm;
    CounterPtr _wLight,_wCharm;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(HRS_1990_I280958);


}
