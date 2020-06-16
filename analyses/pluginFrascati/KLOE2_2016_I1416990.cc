// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta -> 3 pi analysis
  class KLOE2_2016_I1416990 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(KLOE2_2016_I1416990);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_dalitz, 1,1,1);
      double step(0.1), x(-0.9);
      for(unsigned int ix=1;ix<18;++ix) {
	Histo1DPtr temp;
	book(temp,2,1,ix);
	_h_dalitz.add( x, x+step, temp);
	x += step;
      }
      book(_norm,"TMP/norm");
    }
    
    void findDecayProducts(const Particle & mother, unsigned int & nstable,
			   Particles &pi0, Particles &pip, Particles &pim) {
      for(const Particle & p : mother.children()) {
        int id = p.pid();
        if (id == PID::PIMINUS ) {
	  pim.push_back(p);
          ++nstable;
	}
        else if (id == PID::PIPLUS) {
	  pip.push_back(p);
          ++nstable;
        }
        else if (id == PID::PI0) {
	  pi0.push_back(p);
          ++nstable;
        }
        else if ( !p.children().empty() ) {
          findDecayProducts(p, nstable, pi0,pip,pim);
        }
        else
          ++nstable;
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over eta mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
	if(p.pid()!=221) continue;
      	Particles pi0, pip, pim;
      	unsigned nstable(0);
      	findDecayProducts(p,nstable,pi0,pip,pim);
      	if(nstable==3 && pi0.size()==1 && pip.size()==1 && pim.size()==1) {
	  // massesx
      	  double meta = p.mass();
      	  double mpip = pip[0].mass();
      	  double mpim = pim[0].mass();
      	  double mpi0 = pi0[0].mass();
	  // kinetic energies
      	  double Q  = meta-mpip-mpim-mpi0;
      	  double Ep = 0.5/meta*(sqr(meta)+sqr(mpip)-(p.momentum()-pip[0].momentum()).mass2())-mpip;
      	  double Em = 0.5/meta*(sqr(meta)+sqr(mpim)-(p.momentum()-pim[0].momentum()).mass2())-mpim;
      	  double E0 = 0.5/meta*(sqr(meta)+sqr(mpi0)-(p.momentum()-pi0[0].momentum()).mass2())-mpi0;
	  // X and Y variables
      	  double X = sqrt(3.)/Q*(Ep-Em);
      	  double Y = 3.*E0/Q-1.;
       	  _dalitz ->fill(X,Y);
	  _h_dalitz.fill(Y,X);
       	  if(fabs(X)<0.03225806451612903 && Y>0. && Y<0.1) _norm->fill();
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_dalitz,0.06451612903225806*0.1/ *_norm);
      for (Histo1DPtr histo : _h_dalitz.histos()) scale(histo,0.06451612903225806/ *_norm);

    }

    //@}


    /// @name Histograms
    //@{
    Histo2DPtr _dalitz;
    BinnedHistogram _h_dalitz;
    CounterPtr _norm;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(KLOE2_2016_I1416990);


}
