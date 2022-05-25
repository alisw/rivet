// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief kaon production at low energies
  class PLUTO_1981_I165122 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1981_I165122);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      
      // Book histograms
      book(_c_hadrons , "/TMP/sigma_hadrons");
      book(_c_muons   , "/TMP/sigma_muons");
      book(_c_kaons   , "/TMP/sigma_kaons");
      book(_c_hadronsY, "/TMP/sigma_hadronsY");
      book(_c_muonsY  , "/TMP/sigma_muonsY");
      book(_c_kaonsY  , "/TMP/sigma_kaonsY");
      if (isCompatibleWithSqrtS(9.4)) {
        book(_h_spectrum1, 5, 1, 1);
      }
      else if (isCompatibleWithSqrtS(30.0, 1E-2)) {
        book(_h_spectrum1, 4, 1, 1);
      }
      book(_h_spectrum2, 6, 1, 1);
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles &kaons, Particles& stable) {
      for(const Particle & p: mother.children()) {
        const int id = p.pid();
        if(id==130 || id ==310) {
          kaons.push_back(p);
        }
        if (id==111 or p.children().empty())
          stable.push_back(p);
        else
          findDecayProducts(p, kaons, stable);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // Find the Upsilons among the unstables
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particles upsilons = ufs.particles(Cuts::pid==553);
      // Continuum
      if (upsilons.empty()) {
        MSG_DEBUG("No Upsilons found => continuum event");
	// final state particles
	const FinalState& fs = apply<FinalState>(event, "FS");
	map<long,int> nCount;
	int ntotal(0);
	for (const Particle& p : fs.particles()) {
	  nCount[p.pid()] += 1;
	  ++ntotal;
	}
	// mu+mu- + photons
	if(nCount[-13]==1 and nCount[13]==1 &&
	   ntotal==2+nCount[22])
	  _c_muons->fill();
	// everything else
	else
	  _c_hadrons->fill();
	// unstable particles
	for (const Particle& p : ufs.particles(Cuts::pid==130 or Cuts::pid==310)) {
	  if(_h_spectrum1) {
	    double xp = p.p3().mod()/meanBeamMom;
	    _h_spectrum1->fill(xp);
	  }
	  _c_kaons->fill();
	}
      }
      else {
        MSG_DEBUG("Upsilons found => resonance event");
        for (const Particle& ups : upsilons) {
          Particles kaons,stable;
          // Find the decay products we want
          findDecayProducts(ups, kaons, stable);
	  // boost to rest frame (if required)
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();

	  map<long,int> nCount;
	  int ntotal(0);
	  for (const Particle& p : stable) {
	    nCount[p.pid()] += 1;
	    ++ntotal;
	  }
	  for( const Particle & kaon : kaons) {
            const FourMomentum p2 = cms_boost.transform(kaon.momentum());
            const double xp = 2.*p2.p3().mod()/mass;
	    _h_spectrum2->fill(xp);
	    _c_kaonsY->fill();
	  }
	  // mu+mu- + photons
	  if(nCount[-13]==1 and nCount[13]==1 &&
	     ntotal==2+nCount[22])
	    _c_muonsY->fill();
	  // everything else
	  else
	    _c_hadronsY->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // energy dependent
      for(unsigned int ix=1;ix<4;++ix) {
	CounterPtr denom = (ix==1 || ix==3 ) ? _c_muons : _c_hadrons;
	Scatter1D R = *_c_kaons/ *denom;
	double              rval = R.point(0).x();
	pair<double,double> rerr = R.point(0).xErrs();
	Scatter2D temphisto(refData(ix, 1, 1));
	Scatter2DPtr mult;
	book(mult,ix, 1, 1);
	for (size_t b = 0; b < temphisto.numPoints(); b++) {
	  const double x  = temphisto.point(b).x();
	  pair<double,double> ex = temphisto.point(b).xErrs();
	  if(x==9.458) {
	    if(_c_kaonsY->val()>0.) {
	      Scatter1D R2;
	      if(ix==1 ) {
		R2 = *_c_kaonsY/ *_c_muonsY;
	      }
	      else if(ix==2) {
		R2 = *_c_kaonsY/ *_c_hadronsY;
	      }
	      else if(ix==3) {
		R2 = *_c_kaonsY/ *_c_muonsY;
	      }
	      mult   ->addPoint(x, R2.point(0).x(), ex, R2.point(0).xErrs());
	    }
	    else {
	      mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	    }
	  }
	  else if (denom->numEntries()>0) {
	    pair<double,double> ex2 = ex;
	    if(ex2.first ==0.) ex2. first=0.0001;
	    if(ex2.second==0.) ex2.second=0.0001;
	    
	    if (inRange(sqrtS()/GeV, x-ex2.first, x+ex2.second)) {
	      mult   ->addPoint(x, rval, ex, rerr);
	    }
	    else {
	      mult   ->addPoint(x, 0., ex, make_pair(0.,.0));
	    }
	  }
	}
      }
      // normalize the spectra if required
      if(_h_spectrum1) {
	scale(_h_spectrum1, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
      }
      if(_h_spectrum2) {
	scale(_h_spectrum2, 1./_c_hadronsY->val());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_spectrum1, _h_spectrum2;
    CounterPtr _c_hadrons, _c_muons, _c_kaons;
    CounterPtr _c_hadronsY, _c_muonsY, _c_kaonsY;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(PLUTO_1981_I165122);


}
