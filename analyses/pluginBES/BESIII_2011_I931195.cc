// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  psi(2S) -> gamma chi_c0,2
  class BESIII_2011_I931195 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2011_I931195);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::pid==10441 || Cuts::pid==445), "UFS");
      declare(FinalState(), "FS");
      for(unsigned int ichi=0;ichi<2;++ichi) {
	for(unsigned int imeson=0;imeson<2;++imeson) {
	  book(_h_thy[ichi][imeson][0],"TMP/h_"+toString(ichi+1)+"_"+toString(imeson+1)+"_gamma",50,-1.,1.);
	  book(_h_thy[ichi][imeson][1],"TMP/h_"+toString(ichi+1)+"_"+toString(imeson+1)+"_meson",50,-1.,1.);
	  book(_h_thy[ichi][imeson][2],"TMP/h_"+toString(ichi+1)+"_"+toString(imeson+1)+"_phi"  ,50,0.,2.*M_PI);
	  for(unsigned int iy=0;iy<3;++iy)
	    book(_h_exp[ichi][imeson][iy],5+ichi,1+imeson,1+iy);
	}
      }
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<3;++iy)
	    book(_c[ix][iy],"TMP/c_"+toString(ix+1)+"_"+toString(iy+1));
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
      static const double cos20=0.9396926207859084;
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
      // loop over chi_c states
      Particle chi;
      bool matched = false;
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p :  ufs.particles()) {
       	if(p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
	if(ncount==1) {
	  matched = true;
	  for(auto const & val : nRes) {
	    if(val.first==PID::PHOTON) {
	      if(val.second!=1) {
	      matched = false;
	      break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    chi=p;
	    break;
	  }
	}
      }
      if(!matched) vetoEvent;
      // have chi_c find psi2S 
      if(chi.parents().empty() || chi.children().size()!=2 ||
	 chi.children()[0].pid() != -chi.children()[1].pid()) vetoEvent;
      Particle psi2S = chi.parents()[0];
      if(psi2S.pid()!=100443 || psi2S.children().size()!=2) vetoEvent;
      // then the first photon
      Particle gamma1;
      if(psi2S.children()[0].pid()==PID::PHOTON)
	gamma1 = psi2S.children()[0];
      else if(psi2S.children()[1].pid()==PID::PHOTON)
	gamma1 = psi2S.children()[1];
      else
	vetoEvent;
      // now the decay products of the chi_c
      Particle mPlus,mMinus;
      bool foundMeson=false;
      for(unsigned int ix=0;ix<2;++ix) {
	if(chi.children()[ix].pid()==PID::PIPLUS ||
	   chi.children()[ix].pid()==PID::KPLUS ) {
	  foundMeson=true;
	  mPlus=chi.children()[ix];
	}
	else if(chi.children()[ix].pid()==PID::PIMINUS ||
		chi.children()[ix].pid()==PID::KMINUS ) {
	  mMinus=chi.children()[ix];
	}
      }
      if(!foundMeson) vetoEvent;
      // cut on photon angles
      Vector3 aGamma = gamma1.p3().unit();
      double cGammaCut = abs(axis.dot(aGamma));
      // type chi state
      unsigned int ichi= chi.pid()==445 ? 0 : 1;
      // type of meson
      unsigned int imeson = mPlus.pid()==PID::PIPLUS ? 0 : 1;
      // first angle of gamma1 w.r.t beam
      double cGamma = axis.dot(gamma1.momentum().p3().unit());
      _h_thy[ichi][imeson][0]->fill(cGamma);
      // axis in the chi frame
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(chi.momentum().betaVec());
      Vector3 e1z = gamma1.momentum().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      FourMomentum pMeson = boost1.transform(mPlus.momentum());
      Vector3 axis1 = pMeson.p3().unit();
      double cMeson = e1z.dot(axis1);
      _h_thy[ichi][imeson][1]->fill(cMeson);
      double phi = atan2(e1y.dot(axis1),e1x.dot(axis1))+M_PI;
      _h_thy[ichi][imeson][2]->fill(phi);
      // moments to extract multipoles for chi_c2
      if(ichi==0) {
	double sGamma = sqrt(1.-sqr(cGamma));
	double sMeson = sqrt(1.-sqr(cMeson));
	double a3 = -3./sqrt(2.)*cos(phi)*sqr(sMeson)*2.*sMeson*cMeson*2.*sGamma*cGamma;
	double a4 = sqrt(3.)*(3.*sqr(cMeson)-1.)*cos(phi)*2.*sMeson*cMeson*2.*sGamma*cGamma;
	double a5 = sqrt(3./2.)*(3.*sqr(cMeson)-1.)*sqr(sGamma)*sqr(sMeson)*(2.*sqr(cos(phi))-1.);
	_c[imeson][0]->fill(a3);
	_c[imeson][1]->fill(a4);
	_c[imeson][2]->fill(a5);
	_c[2][0]->fill(a3);
	_c[2][1]->fill(a4);
	_c[2][2]->fill(a5);
      }
      // now fill experimental plots with cuts
      if(cGammaCut>0.92 || (cGammaCut>0.8 && cGammaCut<0.86)) vetoEvent;
      // cut on charged particles
      if(abs(axis.dot(mPlus .p3().unit()))>0.93) vetoEvent;
      if(abs(axis.dot(mMinus.p3().unit()))>0.93) vetoEvent;
      // cut on angle of photon w.r.t. charged particles
      if(abs(aGamma.dot(mPlus .p3().unit()))>cos20) vetoEvent;
      if(abs(aGamma.dot(mMinus.p3().unit()))>cos20) vetoEvent;
      // fill histos
      _h_exp[ichi][imeson][0]->fill(cGamma);
      _h_exp[ichi][imeson][1]->fill(cMeson);
      _h_exp[ichi][imeson][2]->fill(phi);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // first normalize the histograms
      for(unsigned int ichi=0;ichi<2;++ichi) {
	for(unsigned int imeson=0;imeson<2;++imeson) {
	  for(unsigned int iy=0;iy<3;++iy) {
	    normalize(_h_thy[ichi][imeson][iy]);
	    normalize(_h_exp[ichi][imeson][iy]);
	  }
	}
      }
      // extract the x and y values
      double x,y;
      pair<double,double> dx,dy;
      for(unsigned int ix=0;ix<3;++ix) {
	Scatter1D X = *_c[ix][0]/ *_c[ix][2];
	Scatter1D Y = *_c[ix][0]/ *_c[ix][1];
	x  = X.point(0).x();
	y  = Y.point(0).x();
	dx = X.point(0).xErrs();
	dy = Y.point(0).xErrs();
	Scatter2DPtr mult;
	book(mult, 1+ix, 1, 1);
	mult->addPoint(0.5,x,make_pair(0.5,0.5),dx);
	book(mult, 1+ix, 1, 2);
	mult->addPoint(0.5,y,make_pair(0.5,0.5),dy);
      }
      // convert x and y to M1 and E2
      double M1 = (3*sqrt(10) + sqrt(30)*x - 2*sqrt(15)*y)/(3.*(sqrt(2) + sqrt(6)*x + 2*sqrt(3)*y));
      double E2 = (2*(2*sqrt(3) - 5*sqrt(2)*y - 2*sqrt(3)*sqr(y) + 4*x*(-1 + sqrt(6)*y)))/
	((sqrt(6) - 6*y)*(sqrt(2) + sqrt(6)*x + 2*sqrt(3)*y));
      double e1 = (-4*sqrt(1.6666666666666667) + 4*sqrt(10)*y)/sqr(sqrt(2) + sqrt(6)*x + 2*sqrt(3)*y);
      double e2 = (-4*sqrt(3.3333333333333335)*(2 + sqrt(3)*x))/sqr(sqrt(2) + sqrt(6)*x + 2*sqrt(3)*y);
      double e3 = (20*(-sqrt(2) + y*(sqrt(3) + 3*sqrt(2)*y)))/((sqrt(6) - 6*y)*sqr(sqrt(2) + sqrt(6)*x + 2*sqrt(3)*y));
      double e4 = (-10*(sqrt(6) - 3*sqrt(2)*x))/(3.*sqr(sqrt(2) + sqrt(6)*x + 2*sqrt(3)*y));
      pair<double,double> dM1 = make_pair(sqrt(sqr(e1*dx.first)+sqr(e2*dy.first)),sqrt(sqr(e1*dx.first)+sqr(e2*dy.first)));
      pair<double,double> dE2 = make_pair(sqrt(sqr(e3*dx.first)+sqr(e4*dy.first)),sqrt(sqr(e3*dx.first)+sqr(e4*dy.first)));
      Scatter2DPtr mult;
      book(mult, 4, 1, 1);
      mult->addPoint(0.5,M1,make_pair(0.5,0.5),dM1);
      book(mult, 4, 1, 2);
      mult->addPoint(0.5,E2,make_pair(0.5,0.5),dE2);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_exp[2][2][3], _h_thy[2][2][3];
    CounterPtr _c[3][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2011_I931195);

}
