// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c*+ spectrum
  class CLEOII_1996_I416471 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CLEOII_1996_I416471);


    /// @name Analysis methods
    ///@{
    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_x,2,1,1);
      book(_r,3,1,1);
      book(_c_xi,"TMP/c_xi");
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

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idXi = 4324;
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.645));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idXi)) {
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
	int sign = p.pid()/p.abspid();
	if(isDecay(p,{sign*4132,sign*211})) {
	  _r->fill(0.5);
	}
      }
      unsigned int nxi = ufs.particles(Cuts::abspid==4132).size();
      _c_xi->fill(nxi);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x,1,false);
      scale(_r, 1./ *_c_xi);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x,_r;
    CounterPtr _c_xi;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(CLEOII_1996_I416471);

}
