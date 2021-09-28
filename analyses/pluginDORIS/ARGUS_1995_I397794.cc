// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Spectrum for D_s2+
  class ARGUS_1995_I397794 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1995_I397794);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_rate1,1,1,1);
      book(_h_rate2,1,2,1);
      book(_h_x,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idDs2 = 435;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.625));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idDs2)) {
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
	int sign = p.pid()/p.abspid();
	if(p.children().size()!=2) continue;
	if(p.children()[0].pid()==sign*421 &&
	   p.children()[1].pid()==sign*321) {
	  _h_rate1->fill(xp);
	  _h_rate2->fill(xp);
	}
	else if(p.children()[1].pid()==sign*421 &&
		p.children()[0].pid()==sign*321) {
	  _h_rate1->fill(xp);
	  _h_rate2->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x);
      scale(_h_rate1,0.3*crossSection()/sumOfWeights()/picobarn);
      scale(_h_rate2,crossSection()/sumOfWeights()/picobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x,_h_rate1,_h_rate2;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(ARGUS_1995_I397794);

}
