// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief p pbar correlations
  class DELPHI_2000_I531568 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_2000_I531568);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const ChargedFinalState cfs;
      declare(cfs        ,"CFS"   );
      declare(Thrust(cfs),"Thrust");
      // book histos
      book(_h_pMp,"_n_pMp",8,0.,1.);
      book(_h_sum,"_n_sum",8,0.,1.);
    }

    void findPP(unsigned int mode, const Vector3 & axis, const Particles & part,
		unsigned int & pp, double & dy) {
      pp=0;
      dy=1e30;
      map<double,Particle> rapOrdered;
      for(const Particle & p : part) {
      	const double mom = dot(axis, p.momentum().p3());
      	const double energy = p.E();
      	const double rap = 0.5 * std::log((energy + mom) / (energy - mom));
      	if(mode==0) {
      	  if(rap>0.)  rapOrdered[rap] = p;
      	}
      	else {
      	  if(rap<=0.) rapOrdered[rap] = p;
      	}
      }
      map<double,Particle>::const_iterator ii[2]={rapOrdered.end(),rapOrdered.end()};
      // map<double,Particle>::const_iterator it,i,im;
      for(map<double,Particle>::const_iterator it=rapOrdered.begin();it!=rapOrdered.end();++it) {
      	if(it->second.abspid()==PID::PROTON) {
	  if(ii[0]!=rapOrdered.end())
	    ii[1]=it;
	  else
	    ii[0]=it;
	}
      }
      // number of particles between the proton and antiproton
      int rank = std::distance(ii[0],ii[1]);
      if(rank>4) return;
      // protom/antiproton next to each other, distance to nearest meson near them
      if(rank==1) {
	map<double,Particle>::const_iterator im=ii[0];--ii[0];
	map<double,Particle>::const_iterator ip=ii[1];++ii[1];
	if(ii[0]!=rapOrdered.begin()) {
	  pp=1;
	  dy = min(2./3.*abs(im->first-ii[0]->first),dy);
	}
	if(ip!=rapOrdered.end()) {
	  pp=1;
	  dy = min(2./3.*abs(ip->first-ii[1]->first),dy);
	}
      }
      else {
	double ycent = 0.5*(ii[0]->first+ii[1]->first);
	double ymin=1e30;
	map<double,Particle>::const_iterator im=rapOrdered.end();
	map<double,Particle>::const_iterator it=ii[0];++it;
	for(;it!=ii[1];++it) {
	  double test = abs(ycent-it->first);
	  if(test<ymin) {
	    im=it;
	    ymin=test;
	  }
	}
	pp=2;
	dy = min(abs(ii[0]->first-im->first),abs(ii[1]->first-im->first));
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Thrust thrust = apply<Thrust>(event,"Thrust");
      Vector3 axis = thrust.thrustAxis();
      const ChargedFinalState cfs = apply<ChargedFinalState>(event,"CFS");
      unsigned int np[2]={0,0}, npbar[2]={0,0};
      for(const Particle & p : cfs.particles()) {
        const double mom = dot(axis, p.momentum().p3());
        const double energy = p.E();
        const double rap = 0.5 * std::log((energy + mom) / (energy - mom));
	if(p.abspid()==PID::PROTON) {
	  unsigned int irap = rap>0 ? 0 : 1;
	  if(p.pid()>0) ++np   [irap];
	  else          ++npbar[irap];
	}
      }
      if(np[0]==1 && npbar[0]==1) {
	unsigned int pp = 0;
	double dy(1e30);
	findPP(0,axis,cfs.particles(),pp,dy);
	if(pp==1) {
	  _h_sum->fill(dy);
	}
	else if(pp==2) {
	  _h_sum->fill(dy);
	  _h_pMp->fill(dy);
	}
      }
      if(np[1]==1 && npbar[1]==1) {
	unsigned int pp = 0;
	double dy(1e30);
	findPP(1,axis,cfs.particles(),pp,dy);
	if(pp==1) {
	  _h_sum->fill(dy);
	}
	else if(pp==2) {
	  _h_sum->fill(dy);
	  _h_pMp->fill(dy);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter2DPtr h_r;
      book(h_r,1,1,1);
      divide(_h_pMp,_h_sum,h_r);
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_sum,_h_pMp;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(DELPHI_2000_I531568);

}
