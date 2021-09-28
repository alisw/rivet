// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief anti-deuteron spectrum in upslion decays
  class BABAR_2014_I1286317 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2014_I1286317);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_p[0],6,1,1);
      book(_h_p[1],6,2,1);
      book(_h_p[2],6,3,1);
      book(_h_p[3],6,4,1);
      book(_h_r[0],1,1,1);
      book(_h_r[1],2,1,1);
      book(_h_r[2],3,1,1);
      book(_h_r[3],5,1,1);
      book(_h_r[4],4,1,1);
      book(_w[0],"TMP/w_0");
      book(_w[1],"TMP/w_1");
      book(_w[2],"TMP/w_2");
      book(_w[3],"TMP/w_3");
    }
    
    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& deut) {
      for(const Particle & p: mother.children()) {
	if(p.pid() == _did) {
	  deut.push_back(p);
	}
	else if(!p.children().empty())
	  findDecayProducts(p, deut);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find upsilon states
      UnstableParticles  ufs = apply<UnstableParticles>(event, "UFS");
      Particles ups = ufs.particles(Cuts::pid==553||Cuts::pid==100553||Cuts::pid==200553);
      // none, then continuum event
      if(ups.empty()) {
	Particles deut = ufs.particles(Cuts::pid==_did);
	_w[3]->fill();
	for(const Particle& p : deut) {
	  double mom = p.momentum().p3().mod();
	  _h_p[3]->fill(mom);
	  _h_r[3]->fill(10.58);
	  _h_r[4]->fill(10.58);
	}
      }
      // upsilon decays
      else {
	for(const Particle & Y : ups ) {
	  unsigned int ihist=2;
	  if(Y.pid()==100553) {
	    ihist=1;
	  }
	  else if(Y.pid()==200553) {
	    ihist=0;
	  }
	  _h_r[ihist]->fill(10.58);
	  _w  [ihist]->fill();
	  Particles deut;
	  findDecayProducts(Y, deut);
	  if(deut.empty()) continue;
	  LorentzTransform boost;
	  if (Y.p3().mod() > 1*MeV)
	    boost = LorentzTransform::mkFrameTransformFromBeta(Y.momentum().betaVec());
	  for(const Particle& p : deut) {
	    double mom = boost.transform(p.momentum()).p3().mod();
	    _h_p[ihist]->fill(mom);
	    _w  [ihist]->fill();
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // upsilon decays
      for(unsigned int ix=0;ix<3;++ix) {
	if(_w[ix]->effNumEntries()<=0.) continue;
	scale(_h_p[ix],1e6/ *_w[ix]);
	scale(_h_r[ix],1./ *_w[ix]);
	Histo1DPtr _h_p[4],_h_r[5];
	CounterPtr _w[4];
      }
      // continuum
      if(_w[3]->effNumEntries()>0.) {
	scale(_h_p[3], crossSection()/sumOfWeights()/femtobarn);
	scale(_h_r[4], crossSection()/sumOfWeights()/femtobarn);	
	scale(_h_r[3],1./ *_w[3]);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_p[4],_h_r[5];
    CounterPtr _w[4];
    
    // deuteron id code
    static const int _did = -1000010020;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2014_I1286317);

}
