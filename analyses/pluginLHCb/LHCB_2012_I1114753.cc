// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi+pi- mass in lambdab* decays
  class LHCB_2012_I1114753 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2012_I1114753);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(),"UFS");
      book(_h_mpipi,1,1,1);
    }

    void findDecayProducts(const Particle & mother,
			   unsigned int & nstable,
			   Particles& pip, Particles& pim,
			   Particles & lambda) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
      	if ( id == PID::PIMINUS) {
	  pim.push_back(p);
	  ++nstable;
	}
       	else if (id == PID::PIPLUS) {
       	  pip.push_back(p);
       	  ++nstable;
       	}
	else if (abs(id)==5122) {
	  lambda.push_back(p);
	  ++nstable;
	}
	else if ( !p.children().empty() ) {
	  findDecayProducts(p,nstable,pip,pim,lambda);
	}
	else
	  ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over unstable particles
      for(const Particle& lb : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==101254)) {
	unsigned int nstable(0);
	Particles pip, pim, lambda;
	findDecayProducts(lb,nstable,pip,pim,lambda);
	// check for lambda
	if(lambda.size() !=1 || nstable !=3 ||
	   pip.size()!=1 || pim.size() !=1 ) continue;
	FourMomentum q = pip[0].momentum()+pim[0].momentum();
	_h_mpipi->fill(q.mass()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mpipi);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mpipi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2012_I1114753);

}
