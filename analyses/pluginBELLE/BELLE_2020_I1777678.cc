// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"

namespace Rivet {


  /// @brief Single and di-hadron spectra
  class BELLE_2020_I1777678 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2020_I1777678);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(ChargedFinalState(Cuts::abspid==211 or Cuts::abspid==321 or Cuts::abspid==2212),"CFS");
      // projections
      FinalState fs;
      declare(fs,"FS");
      declare(Thrust(fs),"Thrust");
      // single particle hists
      vector<int> pdg={211,321,2212};
      for(unsigned int ix=0;ix<3;++ix) {
	book(_s_all   [pdg[ix]],1,ix+1,1);
	book(_s_strong[pdg[ix]],1,ix+1,2);
      }
      // dihadron histograms
      double bins[17]={0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,
      		       0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.00};
      unsigned int i0=1;
      for(unsigned int defn=0;defn<3;++defn) {
	for(unsigned int hemi=0;hemi<3;++hemi) {
	  for(unsigned int ip=0;ip<6;++ip) {
	    i0+=1;
	    unsigned int ymax=16;
	    if(i0==7 || i0==19)      ymax=15;
	    else if(i0>=8 && i0<=12) ymax=14;
	    else if(i0==13)          ymax=13;
	    else if(i0==26)          ymax=10;
	    else if(i0==27||i0==30)  ymax= 9;
	    else if(i0==28)          ymax= 7;
	    else if(i0==29)          ymax= 8;
	    else if(i0==31||i0==44)  ymax= 6;
	    else if(i0==45||i0==48)  ymax= 5;
	    else if(i0==46||i0==47)  ymax= 4;
	    else if(i0==49        )  ymax= 3;
	    for(unsigned int iy=0;iy<ymax;++iy) {
	      Histo1DPtr temp;
	      book(temp,i0,1,iy+1);
	      _d_all   [ip][defn][hemi].add(bins[iy],bins[iy+1],temp);
	      book(temp,i0,2,iy+1);
	      _d_strong[ip][defn][hemi].add(bins[iy],bins[iy+1],temp);
	    }
	  }
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
    
    void fillHistos(int ip,bool strong,bool same,bool opp,
		    const Particle & p1, const Particle & p2) {
      for(unsigned int def=0;def<3;++def) {
	double z1 = 0., z2 = 0.;
	if(def==0) {
	  z1 = 2.*p1.momentum().t()/sqrtS();
	  z2 = 2.*p2.momentum().t()/sqrtS();
	}
	else if(def==1) {
	  z1 = 2.*p1.momentum().t()/sqrtS();
	  z2 = (p1.momentum()*p2.momentum())/p1.momentum().t()/sqrtS();
	}
	else if(def==2) {
	  double p1p2 = p1.momentum()*p2.momentum();
	  double p1q = p1.momentum().t()*sqrtS();
	  double p2q = p2.momentum().t()*sqrtS();
	  z1 = (p1p2-p1.mass2()*p2.mass2()/p1p2)/(p2q-p2.mass2()*p1q/p1p2);
	  z2 = (p1p2-p1.mass2()*p2.mass2()/p1p2)/(p1q-p1.mass2()*p2q/p1p2);
	}
	_d_all[ip][def][0].fill(z1,z2,0.5);
	if(strong) _d_strong[ip][def][0].fill(z1,z2,0.5);
	if(same) {
	  _d_all[ip][def][1].fill(z1,z2,0.5);
	  if(strong) _d_strong[ip][def][1].fill(z1,z2,0.5);
	}
	if(opp) {
	  _d_all[ip][def][2].fill(z1,z2,0.5);
	  if(strong) _d_strong[ip][def][2].fill(z1,z2,0.5);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // apply projection
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      // fill single particle histos
      for (const Particle& p : cfs.particles()) {
	const double z = 2.*p.momentum().t()/sqrtS();
	_s_all   [p.abspid()]->fill(z);
	if(!isWeak(p)) _s_strong[p.abspid()]->fill(z);
      }
      // get thrust
      const Thrust thrust = apply<Thrust>(event,"Thrust");
      ThreeVector axis = thrust.thrustAxis();
      Particles piK = cfs.particles(Cuts::abspid==PID::KPLUS or
				    Cuts::abspid==PID::PIPLUS);
      for(unsigned int ix=0;ix<piK.size();++ix) {
	double dot1 = axis.dot(piK[ix].momentum().p3());
	bool weak1 = isWeak(piK[ix]);
	for(unsigned int iy=0;iy<piK.size();++iy) {
	  if(ix==iy) continue;
	  double dot2 = axis.dot(piK[iy].momentum().p3());
	  bool weak2 = isWeak(piK[iy]);
	  bool strong = !weak1 && !weak2;
	  bool same = thrust.thrust()>0.8 && dot1*dot2>0.;
	  bool opp  = thrust.thrust()>0.8 && dot1*dot2<0.;
	  unsigned int ip=0;
	  if(piK[ix].pid()==PID::PIPLUS) {
	    if(piK[iy].pid()==PID::PIPLUS)       ip=1;
	    else if(piK[iy].pid()==PID::PIMINUS) ip=0;
	    else if(piK[iy].pid()==PID::KPLUS  ) ip=3;
	    else if(piK[iy].pid()==PID::KMINUS)  ip=2;
	  }
	  else if(piK[ix].pid()==PID::PIMINUS) {
	    if(piK[iy].pid()==PID::PIPLUS)       ip=0;
	    else if(piK[iy].pid()==PID::PIMINUS) ip=1;
	    else if(piK[iy].pid()==PID::KPLUS  ) ip=2;
	    else if(piK[iy].pid()==PID::KMINUS ) ip=3;
	  }
	  else if(piK[ix].pid()==PID::KPLUS) {
	    if(piK[iy].pid()==PID::PIPLUS)       ip=3;
	    else if(piK[iy].pid()==PID::PIMINUS) ip=2;
	    else if(piK[iy].pid()==PID::KPLUS)   ip=5;
	    else if(piK[iy].pid()==PID::KMINUS)  ip=4;
	  }
	  else if(piK[ix].pid()==PID::KMINUS) {
	    if(piK[iy].pid()==PID::PIPLUS)       ip=2;
	    else if(piK[iy].pid()==PID::PIMINUS) ip=3;
	    else if(piK[iy].pid()==PID::KPLUS)   ip=4;
	    else if(piK[iy].pid()==PID::KMINUS)  ip=5;
	  }
	  fillHistos(ip,strong,same,opp,piK[ix],piK[iy]);
	}
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for (const auto& kv : _s_all)
	scale(kv.second,crossSection()/femtobarn/sumOfWeights());
      for (const auto& kv : _s_strong)
	scale(kv.second,crossSection()/femtobarn/sumOfWeights());
      for(unsigned int ix=0;ix<6;++ix) {
	for(unsigned int iy=0;iy<3;++iy) {
	  for(unsigned int iz=0;iz<3;++iz) {
	    _d_all   [ix][iy][iz].scale(crossSection()/femtobarn/sumOfWeights(),this);
	    _d_strong[ix][iy][iz].scale(crossSection()/femtobarn/sumOfWeights(),this);
	  }
	}
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    map<int,Histo1DPtr> _s_all,_s_strong;
    BinnedHistogram _d_all[6][3][3],_d_strong[6][3][3];
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BELLE_2020_I1777678);

}
