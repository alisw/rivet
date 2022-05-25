// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c+ spectrum
  class ARGUS_1990_I296522 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1990_I296522);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_obs1,1,1,1);
      book(_h_obs2,1,1,2);
      book(_h_obs3,1,1,3);
      book(_h_all1,1,2,1);
      book(_h_all2,1,2,2);
      book(_h_all3,1,2,3);
      book(_h_x,2,1,1);
    }

    void findDecayProducts(Particle parent, Particles & Xi, Particles & pions,unsigned int & nstable) {
      for(const Particle & p : parent.children()) {
	if(p.abspid()==PID::XIMINUS) {
	  Xi.push_back(p);
	  ++nstable;
	}
	else if(p.abspid()==PID::PIPLUS) {
	  pions.push_back(p);
	  ++nstable;
	}
	else if(!p.children().empty())
	  findDecayProducts(p,Xi,pions,nstable);
	else
	  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const int idXip = 4232,idXi0 = 4132;
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      const double Pmax = sqrt(sqr(Emax)-sqr(2.468));
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==idXip or Cuts::abspid==idXi0)) {
	double xp = p.momentum().p3().mod()/Pmax;
	_h_x->fill(xp);
	Particles Xi,pions;
	unsigned int nstable(0);
	findDecayProducts(p,Xi,pions,nstable);
	if(nstable==2&&Xi.size()==1&&pions.size()==1) {
	  _h_obs1->fill(xp);
	  _h_all1->fill(xp);
	}
	else if(nstable==3&&Xi.size()==1&&pions.size()==2) {
	  _h_obs3->fill(xp);
	  _h_all3->fill(xp);
	}
	else if(nstable==4&&Xi.size()==1&&pions.size()==3) {
	  _h_obs2->fill(xp);
	  _h_all2->fill(xp);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_x);
      scale(_h_obs1,0.5*crossSection()/picobarn/sumOfWeights());
      scale(_h_obs2,0.5*crossSection()/picobarn/sumOfWeights());
      scale(_h_obs3,0.5*crossSection()/picobarn/sumOfWeights());
      scale(_h_all1,    crossSection()/picobarn/sumOfWeights());
      scale(_h_all2,    crossSection()/picobarn/sumOfWeights());
      scale(_h_all3,    crossSection()/picobarn/sumOfWeights());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_obs1,_h_obs2,_h_obs3,_h_all1,_h_all2,_h_all3;
    Histo1DPtr _h_x;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1990_I296522);

}
