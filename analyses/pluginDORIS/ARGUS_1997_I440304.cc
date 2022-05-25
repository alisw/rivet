// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Excited Lambda_c spectra
  class ARGUS_1997_I440304 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1997_I440304);


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
      static const int id2595 = 14122;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.595));
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==id2595)) {
	// spectrum
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
	Particles Lambda_c,pions;
	unsigned int nstable(0);
	findDecayProducts(p,Lambda_c,pions,nstable);
	if(nstable==3&&pions.size()==2&&Lambda_c.size()==1) {
	  _h_rate1->fill(xp);
	  _h_rate2->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x);
      // br for lambda_c mode from pdg 2018
      double br = 0.0623;
      scale(_h_rate1,0.3*br*crossSection()/sumOfWeights()/picobarn);
      scale(_h_rate2,    br*crossSection()/sumOfWeights()/picobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_x,_h_rate1,_h_rate2;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1997_I440304);

}
