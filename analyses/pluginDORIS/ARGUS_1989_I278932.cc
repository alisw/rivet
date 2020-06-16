// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Delta++ rate and spectrum
  class ARGUS_1989_I278932 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1989_I278932);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book( _r_ups,1,1,2);
      book(_r_cont,2,1,2);
      book(_h_p,3,1,1);
      book(_w_ups ,"TMP/w_ups" );
      book(_w_cont,"TMP/w_cont");
    }


    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& delta) {
      // dleta pdg code
      static const long id = 2224;
      for(const Particle & p: mother.children()) {
	if(p.abspid() == id) {
	  delta.push_back(p);
	}
	else if(!p.children().empty())
	  findDecayProducts(p, delta);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      UnstableParticles ufs = apply<UnstableParticles>(event, "UFS");
      Particles ups = ufs.particles(Cuts::pid==553);
      // continuum
      if(ups.empty()) {
	_w_cont->fill();
	_r_cont->fill(sqrtS(),ufs.particles(Cuts::abspid==2224).size());
      }
      // upsilon decays
      else {
	for(const Particle & p : ups) {
	  _w_ups->fill();
	  Particles delta;
	  findDecayProducts(p,delta);
	  if(delta.empty()) continue;
	  LorentzTransform boost;
	  if (p.p3().mod() > 1*MeV)
	    boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	  for(const Particle& del : delta) {
	    _r_ups->fill(9.46);
	    double mom = boost.transform(del.momentum()).p3().mod();
	    _h_p->fill(mom);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_w_cont->effNumEntries()>0) {
	scale(_r_cont, 1./ *_w_cont);
      }
      // direct upsilon decays, i.e. to gg gamma and ggg so need to renormalize
      // factor and values from arxiv:0612019
      if(_w_ups->effNumEntries()>0) {
	double Bmumu=0.0249, rhad = 3.56;
	double fact = 1.-(3.+rhad)*Bmumu;
	scale(_r_ups, 1./fact / *_w_ups);
	scale(_h_p  , 100./fact / *_w_ups);
      }	
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _r_ups,_r_cont;
    Histo1DPtr _h_p;
    CounterPtr _w_ups,_w_cont;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(ARGUS_1989_I278932);

}
