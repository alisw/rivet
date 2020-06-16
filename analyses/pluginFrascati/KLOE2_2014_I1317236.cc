// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  form factor for phi-> eta
  class KLOE2_2014_I1317236 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(KLOE2_2014_I1317236);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      // Book histograms
      book(_h_m, 1, 1, 1);
      book(_weight,"TMP/weight");
    }

    void findDecayProducts(const Particle & mother, unsigned int & nstable, unsigned int & neta, 
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
        else if (id == PID::ETA) {
	  ++neta;
          ++nstable;
        }
	else if (id == PID::GAMMA) {
	  ++ngamma;
	  ++nstable;
	}
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, neta,nep,nem,ngamma,ptot);
        }
        else
          ++nstable;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static double me   = 0.5109989461*MeV;
      static double mphi = 1019.461*MeV;
      static double meta  = 547.862*MeV;

      // Loop over phi mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==333)) {
	unsigned nstable(0),neta(0),nep(0),nem(0),ngamma(0);
	FourMomentum ptot;
	findDecayProducts(p,nstable,neta,nep,nem,ngamma,ptot);
	if(nstable==3 && nem==1 && nem==1 && neta==1) {
	  double q = ptot.mass();
	  double beta = sqrt(1.-4*sqr(me/q));
	  double p = sqrt(sqr(1.+sqr(q)/(sqr(mphi)-sqr(meta)))-4.*sqr(mphi*q/(sqr(mphi)-sqr(meta))));
	  double fact = beta*MeV/q*(1.+2.*sqr(me/q))*pow(p,3);
	  _h_m->fill(q/MeV,1./fact);
	}
	else if(nstable==2 && ngamma ==1 && neta==1) {
	  _weight->fill();
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      static double alpha= 7.2973525664e-3;
      scale(_h_m, 1.5*M_PI/alpha/ *_weight);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_m;
    CounterPtr _weight;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(KLOE2_2014_I1317236);


}
