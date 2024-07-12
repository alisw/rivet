// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief baryon rates and spectra
  class BELLE_2017_I1606201 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2017_I1606201);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      for(unsigned int ix=1;ix<16;++ix) {
	book(_h[ix], ix, 1,  1);
	book(_r[ix], 16, 1, ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles()) {
	const Vector3 mom3 = p.p3();
	double pp = mom3.mod();
	double xp = pp/sqrt(0.25*sqr(sqrtS())-sqr(p.mass()));
	int id = abs(p.pid());
	if(id==3122) {
	  _h[ 1]->fill(xp);
	  _r[ 1]->fill(0.5);
	}
	else if(id==3212) {
	  _h[ 2]->fill(xp);
	  _r[ 2]->fill(0.5);
	}
	else if(id==3224) {
	  _h[ 3]->fill(xp);
	  _r[ 3]->fill(0.5);
	}
	else if(id==3124) {
	  _h[ 4]->fill(xp);
	  _r[ 4]->fill(0.5);
	}
	else if(id==3312) {
	  _h[ 5]->fill(xp);
	  _r[ 5]->fill(0.5);
	}
	else if(id==3334) {
	  _h[ 6]->fill(xp);
	  _r[ 6]->fill(0.5);
	}
	else if(id==3324) {
	  _h[ 7]->fill(xp);
	  _r[ 7]->fill(0.5);
	}
	else if(id==4122) {
	  _h[ 8]->fill(xp);
	  _r[ 8]->fill(0.5);
	}
	else if(id==14122) {
	  _h[ 9]->fill(xp);
	  _r[ 9]->fill(0.5);
	}
	else if(id==4124) {
	  _h[10]->fill(xp);
	  _r[10]->fill(0.5);
	}
	else if(id==4112) {
	  _h[11]->fill(xp);
	  _r[11]->fill(0.5);
	}
	else if(id==4114) {
	  _h[12]->fill(xp);
	  _r[12]->fill(0.5);
	}
	else if(id==4332) {
	  if(isDecay(p,{3334,-211})) {
	    _h[13]->fill(xp);
	    _r[13]->fill(0.5);
	  }
	}
	else if(id==4132) {
	  if(isDecay(p,{3312,211})) {
	    _h[14]->fill(xp);
	    _r[14]->fill(0.5);
	  }
	  else if(isDecay(p,{3334,321})) {
	    _h[15]->fill(xp);
	    _r[15]->fill(0.5);
	  }
	}
      }
    }

    // Check for explicit decay into pdgids
    bool isDecay(const Particle& mother, vector<int> ids) {
      if(mother.pid()<0) {
	for(unsigned int ix=0;ix<ids.size();++ix)
	  ids[ix] *= -1;
      }
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // norm to cross section
      for(unsigned int ix=1;ix<16;++ix) {
	if( ix<=4 || (ix>=8 &&ix<=12) ) 
	  scale(_h[ix], crossSection()/nanobarn/sumOfWeights());
	else
	  scale(_h[ix], crossSection()/picobarn/sumOfWeights());
	scale(_r[ix], crossSection()/picobarn/sumOfWeights());
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h[16],_r[16];
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BELLE_2017_I1606201);


}
