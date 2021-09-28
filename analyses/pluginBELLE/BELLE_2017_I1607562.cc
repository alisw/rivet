// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief BELLE double differential cross section
  class BELLE_2017_I1607562 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2017_I1607562);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      FinalState fs;
      declare(fs,"FS");
      declare(Thrust(fs),"Thrust");
      // histograms
      double bins[17]={0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,
		       0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
      for(unsigned int ip=0;ip<6;++ip) {
	for(unsigned int iy=0;iy<16;++iy) {
	  Histo1DPtr temp;
	  book(temp,1,ip+1,iy+1);
	  _h_all   [ip].add(bins[iy],bins[iy+1],temp);
	  book(temp,2,ip+1,iy+1);
	  _h_strong[ip].add(bins[iy],bins[iy+1],temp);
	}
      }
    }

    bool isWeak(const Particle & p) {
      bool weak = false;
      if(p.parents().empty()) return weak;
      Particle parent = p.parents()[0];
      while (!parent.parents().empty()) {
	if(parent.abspid()==411  || parent.abspid()==421  || parent.abspid()==431  ||
	   parent.abspid()==4122 || parent.abspid()==4232 || parent.abspid()==4132 ||
	   parent.abspid()==4332) {
	  weak=true;
	  break;
	}
	parent = parent.parents()[0];
      }
      return weak;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get thrust and apply cut
      const Thrust thrust = apply<Thrust>(event,"Thrust");
      if(thrust.thrust()<0.8) vetoEvent;
      // get thrust axis
      Vector3 axis = thrust.thrustAxis();
      Particles charged = apply<FinalState>(event,"FS").particles(Cuts::abspid==PID::KPLUS or
								  Cuts::abspid==PID::PIPLUS);
      for(unsigned int ix=0;ix<charged.size();++ix) {
	double dot1 = axis.dot(charged[ix].momentum().p3());
	bool weak1 = isWeak(charged[ix]);
	if(2.*charged[ix].momentum().t()/sqrtS()<0.1) continue;
	for(unsigned int iy=ix+1;iy<charged.size();++iy) {
	  if(2.*charged[iy].momentum().t()/sqrtS()<0.1) continue;
	  double dot2 = axis.dot(charged[iy].momentum().p3());
	  bool weak2 = isWeak(charged[iy]);
	  if(dot1*dot2<0.) continue;
	  FourMomentum p = charged[ix].momentum()+charged[iy].momentum();
	  double z12 = 2.*p.t()/sqrtS();
	  double m12 = p.mass();
	  bool strong = !weak1 && !weak2;
	  if(charged[ix].pid()==PID::PIPLUS) {
	    if(charged[iy].pid()==PID::PIPLUS) {
	      _h_all[1].fill(z12,m12);
	      if(strong) _h_strong[1].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::PIMINUS) {
	      _h_all[0].fill(z12,m12);
	      if(strong) _h_strong[0].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KPLUS) {
	      _h_all[3].fill(z12,m12);
	      if(strong) _h_strong[3].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KMINUS) {
	      _h_all[2].fill(z12,m12);
	      if(strong) _h_strong[2].fill(z12,m12);
	    }
	  }
	  else if(charged[ix].pid()==PID::PIMINUS) {
	    if(charged[iy].pid()==PID::PIPLUS) {
	      _h_all[0].fill(z12,m12);
	      if(strong) _h_strong[0].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::PIMINUS) {
	      _h_all[1].fill(z12,m12);
	      if(strong) _h_strong[1].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KPLUS) {
	      _h_all[2].fill(z12,m12);
	      if(strong) _h_strong[2].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KMINUS) {
	      _h_all[3].fill(z12,m12);
	      if(strong) _h_strong[3].fill(z12,m12);
	    }
	  }
	  else if(charged[ix].pid()==PID::KPLUS) {
	    if(charged[iy].pid()==PID::PIPLUS) {
	      _h_all[3].fill(z12,m12);
	      if(strong) _h_strong[3].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::PIMINUS) {
	      _h_all[2].fill(z12,m12);
	      if(strong) _h_strong[2].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KPLUS) {
	      _h_all[5].fill(z12,m12);
	      if(strong) _h_strong[5].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KMINUS) {
	      _h_all[4].fill(z12,m12);
	      if(strong) _h_strong[4].fill(z12,m12);
	    }
	  }
	  else if(charged[ix].pid()==PID::KMINUS) {
	    if(charged[iy].pid()==PID::PIPLUS) {
	      _h_all[2].fill(z12,m12);
	      if(strong) _h_strong[2].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::PIMINUS) {
	      _h_all[3].fill(z12,m12);
	      if(strong) _h_strong[3].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KPLUS) {
	      _h_all[4].fill(z12,m12);
	      if(strong) _h_strong[4].fill(z12,m12);
	    }
	    else if(charged[iy].pid()==PID::KMINUS) {
	      _h_all[5].fill(z12,m12);
	      if(strong) _h_strong[5].fill(z12,m12);
	    }
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<6;++ix) {
	_h_all   [ix].scale(fact,this);
	_h_strong[ix].scale(fact,this);
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistogram _h_all[6],_h_strong[6];
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2017_I1607562);

}
