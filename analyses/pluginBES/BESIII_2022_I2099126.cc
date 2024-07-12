// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief JPsi > Lambda, Lambdabar with Lambda -> n gamma
  class BESIII_2022_I2099126 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2099126);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");
      for(unsigned int ix=0;ix<2;++ix) {
	book(_n[ix],"TMP/n_" + toString(ix+1));
	book(_t[ix],"TMP/t_" + toString(ix+1));
	for(unsigned int iy=0;iy<2;++iy) {
	  book(_h_mu[ix][iy],1,1,2*ix+iy+1);
	}
      }
      book(_n[2],"TMP/n_3");
      book(_t[2],"TMP/t_3");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // loop over lambda0 baryons
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      Particle Lambda,LamBar;
      bool matched(false);
      for (const Particle& p :  ufs.particles(Cuts::abspid==3122)) {
       	if(p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
       	matched=false;
       	// check for antiparticle
      	for (const Particle& p2 :  ufs.particles(Cuts::pid==-p.pid())) {
      	  if(p2.children().empty()) continue;
      	  map<long,int> nRes2=nRes;
      	  int ncount2 = ncount;
      	  findChildren(p2,nRes2,ncount2);
      	  if(ncount2==0) {
      	    matched = true;
      	    for(auto const & val : nRes2) {
      	      if(val.second!=0) {
      		matched = false;
      		break;
      	      }
      	    }
      	    // found baryon and antibaryon
      	    if(matched) {
	      if(p.pid()>0) {
		Lambda = p;
		LamBar = p2;
	      }
	      else {
		Lambda = p2;
		LamBar = p;
	      }	
       	      break;
       	    }
       	  }
       	}
      	if(matched) break;
      }
      if(!matched) vetoEvent;
      // check the Lambda decay mode
      bool radiative[2]={false,false};
      // identifyt Lambda decay
      Particle baryon1;
      if ( (Lambda.children()[0].pid()==PID::PROTON &&
	    Lambda.children()[1].pid()==PID::PIMINUS ) ) {
	radiative[0]=false;
	baryon1 = Lambda.children()[0];
      }
      else if ( (Lambda.children()[1].pid()==PID::PROTON &&
		 Lambda.children()[0].pid()==PID::PIMINUS ) ) {
	radiative[0]=false;
	baryon1 = Lambda.children()[1];
      }
      else if ( (Lambda.children()[0].pid()==PID::NEUTRON &&
		 Lambda.children()[1].pid()==PID::PHOTON ) ) {
	radiative[0]=true;
	baryon1 = Lambda.children()[0];
      }
      else if ( (Lambda.children()[1].pid()==PID::NEUTRON &&
		 Lambda.children()[0].pid()==PID::PHOTON ) ) {
	radiative[0]=true;
	baryon1 = Lambda.children()[1];
      }
      else
	vetoEvent;
      Particle baryon2;
      if ( (LamBar.children()[0].pid()==PID::ANTIPROTON &&
	    LamBar.children()[1].pid()==PID::PIPLUS ) ) {
	radiative[1]=false;
	baryon2 = LamBar.children()[0];
      }
      else if ( (LamBar.children()[1].pid()==PID::ANTIPROTON &&
		 LamBar.children()[0].pid()==PID::PIPLUS ) ) {
	radiative[1]=false;
	baryon2 = LamBar.children()[1];
      }
      else if ( (LamBar.children()[0].pid()==PID::ANTINEUTRON &&
		 LamBar.children()[1].pid()==PID::PHOTON ) ) {
	radiative[1]=true;
	baryon2 = LamBar.children()[0];
      }
      else if ( (LamBar.children()[1].pid()==PID::ANTINEUTRON &&
		 LamBar.children()[0].pid()==PID::PHOTON ) ) {
	radiative[1]=true;
	baryon2 = LamBar.children()[1];
      }
      else
	vetoEvent;
      if (radiative[0] == radiative[1]) vetoEvent;
      // boost to the Lambda rest frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambda.momentum().betaVec());
      Vector3 e1z = Lambda.momentum().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      Vector3 axis1 = boost1.transform(baryon1.momentum()).p3().unit();
      double n1x(e1x.dot(axis1)),n1y(e1y.dot(axis1)),n1z(e1z.dot(axis1));
      // boost to the Lambda bar
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(LamBar.momentum().betaVec());
      Vector3 axis2 = boost2.transform(baryon2.momentum()).p3().unit();
      double n2x(e1x.dot(axis2)),n2y(e1y.dot(axis2)),n2z(e1z.dot(axis2));
      double cosL = axis.dot(Lambda.momentum().p3().unit());
      double sinL = sqrt(1.-sqr(cosL));
      double T1 = sqr(sinL)*n1x*n2x+sqr(cosL)*n1z*n2z;
      // lambda -> n gamma
      if(radiative[0]) {
	_h_mu[0][0]->fill( cosL,n2y);
	_h_mu[0][1]->fill( cosL,n1y);
	_n[0]->fill();
	_n[2]->fill();
	_t[0]->fill(T1);
	_t[2]->fill(T1);
      }
      // lambdabar -> nbar gamma
      else {
	_h_mu[1][0]->fill( cosL,n1y);
	_h_mu[1][1]->fill( cosL,n2y);
	_n[1]->fill();
	_n[2]->fill();
	_t[1]->fill(T1);
	_t[2]->fill(T1);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // values of constants
      double aPsi  = 0.461;
      double aPlus =-0.758;
      double factor = 45.*(3. +aPsi)/(11. + 5.*aPsi)/aPlus;
      // plots
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  scale(_h_mu[ix][iy],10.*0.2/ *_n[ix]);
	}
      }
      // alpha from the moments
      for(unsigned int ix=0;ix<3;++ix) {
	double value = _t[ix]->val()/_n[ix]->val();
	double error = _t[ix]->err()/_n[ix]->val();
 	value *= factor;
	error *= abs(factor);
	if(ix==1) value *=-1.;
	Scatter2DPtr  alpha;
	book(alpha,2,1,1+ix);
	alpha->addPoint(0.5, value, make_pair(0.5,0.5),
			make_pair(error,error) );
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mu[2][2];
    CounterPtr _n[3];
    Histo1DPtr _h_ctheta[3];
    CounterPtr _t[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2099126);

}
