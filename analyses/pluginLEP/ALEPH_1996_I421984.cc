// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh" 

namespace Rivet {


  /// @brief tau -> omega pi(pi)
  class ALEPH_1996_I421984 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALEPH_1996_I421984);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_omegapi_momega[0], 1, 1, 1);
      book(_h_omegapi_momega[1], 2, 1, 1);
      book(_h_omegapi_momegapi,  3, 1, 1);
      book(_h_omegapipi_momega , 6, 1, 1);
    }

    void findDecayProducts(const Particle &mother, unsigned int & nstable,
                           Particles& pip, Particles& pim, Particles & pi0) {
      for (const Particle &p : mother.children()) {
        long id = p.pid();
        if (id == PID::PI0 ) {
	  pi0.push_back(p);
          ++nstable;
	}
        else if (id == PID::K0S || id == PID::K0L)
          ++nstable;
        else if (id == PID::PIPLUS) {
          pip.push_back(p);
          ++nstable;
        }
        else if (id == PID::PIMINUS) {
          pim.push_back(p);
          ++nstable;
        }
        else if (id == PID::KPLUS) {
          ++nstable;
        }
        else if (id == PID::KMINUS) {
          ++nstable;
        }
        else if (!p.children().empty()) {
          findDecayProducts(p, nstable, pip, pim, pi0);
        }
        else  ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Find the taus
      Particles taus;
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::TAU)) {
        Particles pip, pim, pi0;
        unsigned int nstable = 0;
        // Find the decay products we want
        findDecayProducts(p, nstable, pip, pim, pi0);
        if (p.pid() < 0) {
          swap(pip, pim);
        }
	if(nstable==5&&pip.size()==1&&pim.size()==2&&pi0.size()==1) {
	  for(unsigned int ix=0;ix<2;++ix) {
	    double momega = (pim[ix].momentum()+pip[0].momentum()+pi0[0].momentum()).mass();
	    _h_omegapi_momega[0]->fill(momega);
	    _h_omegapi_momega[1]->fill(momega);
	  }
	  double mtotal = (pim[0].momentum()+pim[1].momentum()+pip[0].momentum()+pi0[0].momentum()).mass();
	  _h_omegapi_momegapi->fill(mtotal);
	}
	else if(nstable==6&&pip.size()==1&&pim.size()==2&&pi0.size()==2) {
	  for(unsigned int ix=0;ix<2;++ix) {
	    for(unsigned int iy=0;iy<2;++iy) {
	      double momega = (pim[ix].momentum()+pip[0].momentum()+pi0[iy].momentum()).mass();
	      _h_omegapipi_momega->fill(momega);
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // normalized measurement
      normalize(_h_omegapi_momegapi);
      // normalize to integral of data
      normalize(_h_omegapi_momega[0],4164.8*0.01 );
      normalize(_h_omegapi_momega[1],5936.3*0.01 );
      normalize(_h_omegapipi_momega ,1967.6*0.035);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_omegapi_momega[2], _h_omegapi_momegapi,_h_omegapipi_momega;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ALEPH_1996_I421984);


}
