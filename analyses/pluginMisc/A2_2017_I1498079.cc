// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief form factors in pi0->gamma gamma(e+e-)
  class A2_2017_I1498079 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(A2_2017_I1498079);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_h_pi0  , 1, 1, 1);
      book(_weight_pi0  ,"TMP/weight_pi0");

    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable,
                           unsigned int & nep, unsigned int & nem, unsigned int & ngamma,
			   FourMomentum & ptot) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::EMINUS ) {
	  ++nem;
          ++nstable;
	  ptot += p.momentum();
	}
        else if (id == PID::EPLUS) {
          ++nep;
          ++nstable;
	  ptot += p.momentum();
        }
	else if (id == PID::GAMMA && p.children().empty() ) {
	  ++ngamma;
	  ++nstable;
	}
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable,nep,nem,ngamma,ptot);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static double me     = 0.5109989461*MeV;
      static double mpi    = 134.9770*MeV;

      // Loop over pi0 and omega mesons
      for( const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==111)) {
	unsigned nstable(0),nep(0),nem(0),ngamma(0);
	FourMomentum ptot;
	findDecayProducts(p,nstable,nep,nem,ngamma,ptot);
	if(nstable==3 && nem==1 && nem==1 && ngamma==1) {
	  double q = ptot.mass();
	  double beta = sqrt(1.-4*sqr(me/q));
	  double p = 1.-sqr(q/mpi);
	  double fact = beta*MeV/q*(1.+2.*sqr(me/q))*pow(p,3);
	  _h_pi0->fill(q/MeV,1./fact);
	}
	else if(nstable==2 && ngamma==2) {
	  _weight_pi0->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      static double alpha= 7.2973525664e-3;
      scale(_h_pi0  , 0.75*M_PI/alpha/ *_weight_pi0);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pi0;
    CounterPtr _weight_pi0;

    ///@}


  };


  RIVET_DECLARE_PLUGIN(A2_2017_I1498079);

}
