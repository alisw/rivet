// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Spectrum for Lambda_c(2625)
  class ARGUS_1993_I357132 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ARGUS_1993_I357132);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_rate,1,1,1);
      book(_h_x,2,1,1);
    }
    
    void findDecayProducts(Particle parent, Particles & Lambda_c, Particles & pions,unsigned int & nstable) {
      for(const Particle & p : parent.children()) {
	if(p.abspid()==4122) {
	  Lambda_c.push_back(p);
	  ++nstable;
	}
	else if(p.abspid()==PID::PIPLUS) {
	  pions.push_back(p);
	  ++nstable;
	}
	else if(!p.children().empty())
	  findDecayProducts(p,Lambda_c,pions,nstable);
	else
	  ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int id2625 = 4124;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.625));
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==id2625)) {
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
	Particles Lambda_c,pions;
	unsigned int nstable(0);
	findDecayProducts(p,Lambda_c,pions,nstable);
	if(nstable==3&&pions.size()==2&&Lambda_c.size()==1)
	  _h_rate->fill(10.58);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x);
      scale(_h_rate,crossSection()/picobarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x,_h_rate;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(ARGUS_1993_I357132);

}
