// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Sigma_c(2800) rate and spectra
  class BELLE_2004_I668024 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I668024);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      for(unsigned int ix=0;ix<3;++ix) {
	book(_r[ix], 1, 1, 1+ix);
	book(_h[ix], 2, 1, 1+ix);
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
	if(id==_sid+110 && isDecay(p,{4122,-211})) {
	  _h[0]->fill(xp);
	  _r[0]->fill(0.5);
	}
	else if(id==_sid+210&& isDecay(p,{4122,111})) {
	  _h[1]->fill(xp);
	  _r[1]->fill(0.5);
	}
	else if(id==_sid+220&& isDecay(p,{4122,211})) {
	  _h[2]->fill(xp);
	  _r[2]->fill(0.5);
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
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h[ix]);
	scale(_r[ix], crossSection()/picobarn/sumOfWeights());
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h[3],_r[3];
    static const int _sid = 14002;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2004_I668024);

}
