// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief form factors for omega->pi and eta->gamma
  class A2_2017_I1486671 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(A2_2017_I1486671);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_eta  , 1, 1, 1);
      book(_h_omega, 2, 1, 1);
      book(_weight_eta  ,"TMP/weight_eta");
      book(_weight_omega,"TMP/weight_omega");
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & npi, 
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
        else if (id == PID::PI0) {
	  ++npi;
          ++nstable;
        }
	else if (id == PID::GAMMA && p.children().empty() ) {
	  ++ngamma;
	  ++nstable;
	}
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, npi,nep,nem,ngamma,ptot);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static double me     = 0.5109989461*MeV;
      static double momega = 1019.461*MeV;
      static double meta   = 547.862 *MeV;
      static double mpi    = 134.9770*MeV;

      // Loop over eta and omega mesons
      for( const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==223 or Cuts::pid==221)) {
       	unsigned nstable(0),npi(0),nep(0),nem(0),ngamma(0);
       	FourMomentum ptot;
      	findDecayProducts(p,nstable,npi,nep,nem,ngamma,ptot);
	if(p.pid()==221) {
	  if(nstable==3 && nem==1 && nem==1 && ngamma==1) {
	    double q = ptot.mass();
	    double beta = sqrt(1.-4*sqr(me/q));
	    double p = 1.-sqr(q/meta);
	    double fact = beta*MeV/q*(1.+2.*sqr(me/q))*pow(p,3);
	    _h_eta->fill(q/MeV,1./fact);
	  }
	  else if(nstable==2 && ngamma==2) {
	    _weight_eta->fill();
	  }
	}
	else {
	  if(nstable==3 && nem==1 && nem==1 && npi==1) {
	    double q = ptot.mass();
	    double beta = sqrt(1.-4*sqr(me/q));
	    double p = sqrt(sqr(1.+sqr(q)/(sqr(momega)-sqr(mpi)))-4.*sqr(momega*q/(sqr(momega)-sqr(mpi))));
	    double fact = beta*MeV/q*(1.+2.*sqr(me/q))*pow(p,3);
	    _h_omega->fill(q/MeV,1./fact);
	  }
	  else if(nstable==2 && ngamma ==1 && npi==1) {
	    _weight_omega->fill();
	  }
	}
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      static double alpha= 7.2973525664e-3;
      scale(_h_eta  , 0.75*M_PI/alpha/ *_weight_eta  );
      scale(_h_omega, 1.5 *M_PI/alpha/ *_weight_omega);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_eta,_h_omega;
    CounterPtr _weight_eta,_weight_omega;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(A2_2017_I1486671);


}
