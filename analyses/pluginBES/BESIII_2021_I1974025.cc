// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > Lambda, Lambdabar
  class BESIII_2021_I1974025 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1974025);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_wsum,"TMP/wsum");
      // for(unsigned int ix=0;ix<6;++ix)
      // 	book(_h_F[ix],1,1,1+ix);
      for(unsigned int ix=0;ix<6;++ix)
      	book(_h_F[ix],"TMP/F_"+toString(ix+1),20,-1.,1.);
      book(_h_F[5],1,1,6);
      book(_h_mu,2,1,1);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // find the Lambdas
      bool matched = false;
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      Particle Lambda,LamBar;
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
	if(abs(p1.pid())!=3122) continue;
	// check fs
	bool fs = true;
	for (const Particle & child : p1.children()) {
	  if(child.pid()==p1.pid()) {
	    fs = false;
	    break;
	  }
	}
	if(!fs) continue;
	// find the children
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
	  matched=false;
	  const Particle& p2 = ufs.particles()[iy];
	  if(abs(p2.pid())!=3122) continue;
	  // check fs
	  bool fs = true;
	  for (const Particle & child : p2.children()) {
	    if(child.pid()==p2.pid()) {
	      fs = false;
	      break;
	    }
	  }
	  if(!fs) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    if(p1.pid()==PID::LAMBDA) {
	      Lambda=p1;
	      LamBar=p2;
	    }
	    else {
	      Lambda=p2;
	      LamBar=p1;
	    }
	    break;
	  }
	}
	if(matched) break;
      }
      // and the children
      Particle proton;
      matched = false;
      for (const Particle & p : Lambda.children()) {
	if(p.pid()==2212) {
	  matched=true;
	  proton=p;
	}
	else if(p.pid()==PID::PHOTON)
	  vetoEvent;
      }
      if(!matched) vetoEvent;
      Particle baryon;
      matched = false;
      for (const Particle & p : LamBar.children()) {
	if(p.pid()==-2212) {
	  baryon=p;
	  matched=true;
	}
	else if(p.pid()==PID::PHOTON)
	  vetoEvent;
      }
      if(!matched) vetoEvent;
      // now for the polarization measurements
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambda.momentum().betaVec());
      Vector3 e1z = Lambda.momentum().p3().unit();
      Vector3 e1y = e1z.cross(axis).unit();
      Vector3 e1x = e1y.cross(e1z).unit();
      Vector3 axis1 = boost1.transform(proton.momentum()).p3().unit();
      double n1x(e1x.dot(axis1)),n1y(e1y.dot(axis1)),n1z(e1z.dot(axis1));
      // boost to the Lambda bar
      LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(LamBar.momentum().betaVec());
      Vector3 axis2 = boost2.transform(baryon.momentum()).p3().unit();
      double n2x(e1x.dot(axis2)),n2y(e1y.dot(axis2)),n2z(e1z.dot(axis2));
      double cosL = -axis.dot(Lambda.momentum().p3().unit());
      double sinL = sqrt(1.-sqr(cosL));
      double T1 = sqr(sinL)*n1x*n2x+sqr(cosL)*n1z*n2z;
      double T2 = sinL*cosL*(n1x*n2z+n1z*n2x);
      double T3 = sinL*cosL*n1y;
      double T4 = sinL*cosL*n2y;
      double T5 = n1z*n2z-sqr(sinL)*n1y*n2y;
      double mu = n1y-n2y;
      _h_F[0]->fill(cosL,T1);
      _h_F[1]->fill(cosL,T2);
      _h_F[2]->fill(cosL,T3);
      _h_F[3]->fill(cosL,T4);
      _h_F[4]->fill(cosL,T5);
      _h_F[5]->fill(cosL);
      _h_mu->fill(cosL,mu);
      _wsum->fill();
    }
    
    pair<double,pair<double,double> > calcAlpha0(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,make_pair(0.,0.));
      double d = 3./(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double c = 3.*(hist->xMax()-hist->xMin())/(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double sum1(0.),sum2(0.),sum3(0.),sum4(0.),sum5(0.);
      for (auto bin : hist->bins() ) {
       	double Oi = bin.area();
	if(Oi==0.) continue;
	double a =  d*(bin.xMax() - bin.xMin());
	double b = d/3.*(pow(bin.xMax(),3) - pow(bin.xMin(),3));
       	double Ei = bin.areaErr();
	sum1 +=   a*Oi/sqr(Ei);
	sum2 +=   b*Oi/sqr(Ei);
	sum3 += sqr(a)/sqr(Ei);
	sum4 += sqr(b)/sqr(Ei);
	sum5 +=    a*b/sqr(Ei);
      }
      // calculate alpha
      double alpha = (-c*sum1 + sqr(c)*sum2 + sum3 - c*sum5)/(sum1 - c*sum2 + c*sum4 - sum5);
      // and error
      double cc = -pow((sum3 + sqr(c)*sum4 - 2*c*sum5),3);
      double bb = -2*sqr(sum3 + sqr(c)*sum4 - 2*c*sum5)*(sum1 - c*sum2 + c*sum4 - sum5);
      double aa =  sqr(sum1 - c*sum2 + c*sum4 - sum5)*(-sum3 - sqr(c)*sum4 + sqr(sum1 - c*sum2 + c*sum4 - sum5) + 2*c*sum5);      
      double dis = sqr(bb)-4.*aa*cc;
      if(dis>0.) {
	dis = sqrt(dis);
	return make_pair(alpha,make_pair(0.5*(-bb+dis)/aa,-0.5*(-bb-dis)/aa));
      }
      else {
	return make_pair(alpha,make_pair(0.,0.));
      }
    }

    pair<double,double> calcCoeff(unsigned int imode,Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai(0.),bi(0.);
	if(imode==0) {
	  bi = (pow(1.-sqr(bin.xMin()),1.5) - pow(1.-sqr(bin.xMax()),1.5))/3.;
	}
	else if(imode>=2 && imode<=4) {
	  bi =   ( pow(bin.xMin(),3)*( -5. + 3.*sqr(bin.xMin()))  +
		   pow(bin.xMax(),3)*(  5. - 3.*sqr(bin.xMax())))/15.;
	}
	else
	  assert(false);
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize histograms
      for(unsigned int ix=0;ix<6;++ix)
	scale(_h_F[ix], 1./ *_wsum);
      scale(_h_mu, 10./ *_wsum);
      // value of aLambda assumed in paper
      double aLambda = 0.754;
      // calculate alpha0
      pair<double,pair<double,double> > alpha0 = calcAlpha0(_h_F[5]);
      Scatter2DPtr _h_alpha0;
      book(_h_alpha0,3,1,1);
      _h_alpha0->addPoint(0.5, alpha0.first, make_pair(0.5,0.5),
			  make_pair(alpha0.second.first,alpha0.second.second) );
      double s2 = -1. + sqr(alpha0.first);
      double s3 = 3 + alpha0.first;
      double s1 = sqr(s3);
      // alpha- and alpha+ from proton data
      pair<double,double> c_T2_p = calcCoeff(2,_h_F[1]); 
      pair<double,double> c_T3_p = calcCoeff(3,_h_F[2]);
      pair<double,double> c_T4_p = calcCoeff(4,_h_F[3]);
      double s4 = sqr(c_T2_p.first);
      double s5 = sqr(c_T3_p.first);
      double s6 = sqr(c_T4_p.first);
      double disc = s1*s5*s6*(-9.*s2*s4 + 4.*s1*s5*s6);
      // now for Delta
      if(disc>0) {
	double sDelta = (-2.*(3. + alpha0.first)*c_T3_p.first)/(aLambda*sqrt(1 - sqr(alpha0.first)));
	double cDelta = (-3*(3 + alpha0.first)*c_T2_p.first)/(-aLambda*aLambda*sqrt(1 - sqr(alpha0.first)));
	double Delta = asin(sDelta);
	if(cDelta<0.) Delta = M_PI-Delta;
	double ds_P = (-9*c_T2_p.first*((-1 + alpha0.first)*(1 + alpha0.first)*  (3 + alpha0.first)*c_T3_p.first*c_T4_p.first*c_T2_p.second +  c_T2_p.first*c_T4_p.first*(c_T3_p.first*(alpha0.second.first + 3*alpha0.first*alpha0.second.first) -(-1 + alpha0.first)*(1 + alpha0.first)*(3 + alpha0.first)*c_T3_p.second)
			      -  (-1 + alpha0.first)*(1 + alpha0.first)*  (3 + alpha0.first)*c_T2_p.first*c_T3_p.first*c_T4_p.second)*disc)/
	  (pow(1 - pow(alpha0.first,2),1.5)*pow(c_T4_p.first,3)*pow(-((disc + 2*s1*s5*s6)/   (s2*s6)),1.5)*(-9*s2*s4 + 4*s1*s5*s6));
	double ds_M = (-9*c_T2_p.first*((-1 + alpha0.first)*(1 + alpha0.first)*  (3 + alpha0.first)*c_T3_p.first*c_T4_p.first*c_T2_p.second +  c_T2_p.first*c_T4_p.first*(c_T3_p.first*(alpha0.second.second + 3*alpha0.first*alpha0.second.second) -(-1 + alpha0.first)*(1 + alpha0.first)*(3 + alpha0.first)*c_T3_p.second)
			      -  (-1 + alpha0.first)*(1 + alpha0.first)*  (3 + alpha0.first)*c_T2_p.first*c_T3_p.first*c_T4_p.second)*disc)/
	  (pow(1 - pow(alpha0.first,2),1.5)*pow(c_T4_p.first,3)*pow(-((disc + 2*s1*s5*s6)/   (s2*s6)),1.5)*(-9*s2*s4 + 4*s1*s5*s6));
	ds_P /= sqrt(1.-sqr(sDelta));
	ds_M /= sqrt(1.-sqr(sDelta));
	Scatter2DPtr _h_sin;
	book(_h_sin,3,1,2);
	_h_sin->addPoint(0.5, Delta/M_PI*180., make_pair(0.5,0.5), make_pair( -ds_P/M_PI*180., -ds_M/M_PI*180. ) );
      }
      // scale to number of observed events in experiment
      // for(unsigned int ix=0;ix<5;++ix)
      // 	scale(_h_F[ix], 262.);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_F[6],_h_mu;
    CounterPtr _wsum;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1974025);

}
