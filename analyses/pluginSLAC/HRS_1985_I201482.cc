// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief event shapes at 29 GeV
  class HRS_1985_I201482 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1985_I201482);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const ChargedFinalState cfs;
      declare(cfs, "FS");
      declare(Sphericity(cfs), "Sphericity");
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");

      // Book histograms
      book(_histSphericity,  1, 1, 1);
      book(_histThrust    ,  3, 1, 1);
      book(_histThrust2Jet,  4, 1, 1);
      book(_histAplanarity,  6, 1, 1);
      book(_histZ         , 10, 1, 1);
      book(_histZ2Jet     , 11, 1, 1);
      book(_histZScale    , 12, 1, 1);
      book(_histZJet[0]   , 13, 1, 1);
      book(_histZJet[1]   , 14, 1, 1);
      book(_histZJet[2]   , 15, 1, 1);
      book(_histXFeyn     , 16, 1, 1);
      book(_histXFeyn2Jet , 17, 1, 1);
      book(_histRap       , 19, 1, 1);
      book(_histRap2Jet   , 20, 1, 1);
      book(_histPtT       , 22, 1, 1);
      book(_histPtT2Jet   , 23, 1, 1);
      book(_histPtTIn     , 24, 1, 1);
      book(_histPtTOut    , 25, 1, 1);
      book(_wSum ,"TMP/wSum");
      book(_wSum2,"TMP/wSum2");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // require 5 charged particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();
      if(numParticles<5) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // calc thrust and sphericity
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      Vector3 axis = thrust.thrustAxis();
      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");
      // identify two and three jet regions
      bool twoJet   = sphericity.sphericity()<=0.25 && sphericity.aplanarity()<=0.1;
      //bool threeJet = sphericity.sphericity() >0.25 && sphericity.aplanarity()<=0.1;
      _wSum->fill();
      if(twoJet) _wSum2->fill();
      // basic event shapes
      _histSphericity->fill(sphericity.sphericity());
      _histThrust    ->fill(thrust.thrust());
      _histAplanarity->fill(sphericity.aplanarity());
      if(twoJet)
	_histThrust2Jet->fill(thrust.thrust());
      double pTSqIn  = 0.;
      double pTSqOut = 0.;
      unsigned int iPlus(0),iMinus(0);
      // single particle  dists
      for(const Particle & p : sortBy(fs.particles(),cmpMomByP)) {
	const double z  = p.p3().mod()/meanBeamMom;
	const double momT = axis.dot(p.p3());
	const double xF = fabs(momT)/meanBeamMom;
        const double energy = p.E();
        const double rap = 0.5 * std::log((energy + momT) / (energy - momT));
        const double pTin  = dot(p.p3(), thrust.thrustMajorAxis());
        const double pTout = dot(p.p3(), thrust.thrustMinorAxis());
	const double pT2 = sqr(pTin)+sqr(pTout);
	pTSqIn  += sqr(dot(p.p3(), sphericity.sphericityMajorAxis()));
	pTSqOut += sqr(dot(p.p3(), sphericity.sphericityMinorAxis()));
	_histZ     ->fill(z         );
	_histZScale->fill(z         );
	_histXFeyn ->fill(xF        ,z);
	_histRap   ->fill(rap       );
	_histPtT   ->fill(pT2       );
	if(twoJet) {
	  _histZ2Jet    ->fill(z  );
	  _histXFeyn2Jet->fill(xF ,z);
	  _histRap2Jet  ->fill(rap);
	  _histPtT2Jet  ->fill(pT2);
	  if(momT>0.&&iPlus<3) {
	    _histZJet[iPlus]->fill(z);
	    iPlus+=1;
	  }
	  else if(momT<0.&&iMinus<3) {
	    _histZJet[iMinus]->fill(z);
	    iMinus+=1;
	  }
	}
      }
      _histPtTIn ->fill(pTSqIn /numParticles);
      _histPtTOut->fill(pTSqOut/numParticles);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_histSphericity);
      normalize(_histThrust);
      normalize(_histThrust2Jet);
      normalize(_histAplanarity);
      scale(_histZ        ,1./ *_wSum);
      scale(_histZScale   , sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      scale(_histXFeyn    ,1./M_PI/ *_wSum);
      scale(_histRap      ,1./ *_wSum);
      scale(_histZ2Jet    ,1./ *_wSum2);
      scale(_histXFeyn2Jet,1./M_PI/ *_wSum2);
      scale(_histRap2Jet  ,1./ *_wSum2);
      scale(_histPtT      ,1./ *_wSum);
      scale(_histPtT2Jet  ,1./ *_wSum2);
      scale(_histPtTIn    ,1./ *_wSum);
      scale(_histPtTOut   ,1./ *_wSum);
      for(unsigned int i=0;i<3;++i)
	scale(_histZJet[i]   ,0.5/ *_wSum2);
    }

    //@}


    /// @name Histograms
    //@{      
    Histo1DPtr _histSphericity, _histThrust, _histThrust2Jet, _histAplanarity,
      _histZ, _histZ2Jet, _histZScale, _histXFeyn, _histXFeyn2Jet, _histRap,
      _histRap2Jet, _histPtT, _histPtT2Jet, _histPtTIn, _histPtTOut ,_histZJet[3];
    CounterPtr _wSum,_wSum2;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(HRS_1985_I201482);


}
