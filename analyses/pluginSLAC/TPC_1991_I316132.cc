// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief D* poliarzation at 29 GeV
  class TPC_1991_I316132 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TPC_1991_I316132);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      {Histo1DPtr tmp; _h_ctheta.add(0.3,0.4,book(tmp, "ctheta_0",20,-1,1));}
      {Histo1DPtr tmp; _h_ctheta.add(0.4,0.5,book(tmp, "ctheta_1",20,-1,1));}
      {Histo1DPtr tmp; _h_ctheta.add(0.5,0.6,book(tmp, "ctheta_2",20,-1,1));}
      {Histo1DPtr tmp; _h_ctheta.add(0.6,0.7,book(tmp, "ctheta_3",20,-1,1));}
      {Histo1DPtr tmp; _h_ctheta.add(0.7,0.8,book(tmp, "ctheta_4",20,-1,1));}
      {Histo1DPtr tmp; _h_ctheta.add(0.8,1.0,book(tmp, "ctheta_5",20,-1,1));}
      book(_h_ctheta_all, "ctheta_all",20,-1,1);

      {Histo1DPtr tmp; _h_phi.add(0.3,0.4,book(tmp, "phi_0",20,-M_PI,M_PI));}
      {Histo1DPtr tmp; _h_phi.add(0.4,0.5,book(tmp, "phi_1",20,-M_PI,M_PI));}
      {Histo1DPtr tmp; _h_phi.add(0.5,0.6,book(tmp, "phi_2",20,-M_PI,M_PI));}
      {Histo1DPtr tmp; _h_phi.add(0.6,0.7,book(tmp, "phi_3",20,-M_PI,M_PI));}
      {Histo1DPtr tmp; _h_phi.add(0.7,0.8,book(tmp, "phi_4",20,-M_PI,M_PI));}
      {Histo1DPtr tmp; _h_phi.add(0.8,1.0,book(tmp, "phi_5",20,-M_PI,M_PI));}
      book(_h_phi_all, "phi_all",20,-M_PI,M_PI);

      {Histo1DPtr tmp; _h_01.add(0.3,0.4,book(tmp, "h_01_0",20,-1,1));}
      {Histo1DPtr tmp; _h_01.add(0.4,0.5,book(tmp, "h_01_1",20,-1,1));}
      {Histo1DPtr tmp; _h_01.add(0.5,0.6,book(tmp, "h_01_2",20,-1,1));}
      {Histo1DPtr tmp; _h_01.add(0.6,0.7,book(tmp, "h_01_3",20,-1,1));}
      {Histo1DPtr tmp; _h_01.add(0.7,0.8,book(tmp, "h_01_4",20,-1,1));}
      {Histo1DPtr tmp; _h_01.add(0.8,1.0,book(tmp, "h_01_5",20,-1,1));}
      book(_h_01_all, "h_01_all",20,-1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double Emax = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();

      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");
      for  (const Particle& p : ufs.particles(Cuts::abspid==413)) {
	if(p.children().size()!=2) continue;
	int sign = p.pid()/413;
	Particle D0;
	if(p.children()[0].pid()==sign*421 && p.children()[1].pid()==sign*211) {
	  D0 = p.children()[0];
	}
	else if(p.children()[1].pid()==sign*421 && p.children()[0].pid()==sign*211) {
	  D0 = p.children()[1];
	}
	else
	  continue;
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	double xE = p.momentum().t()/Emax;
	Vector3 e1z = p.momentum().p3().unit();
	Vector3 e1y = e1z.cross(axis).unit();
	Vector3 e1x = e1y.cross(e1z).unit();
	Vector3 axis1 = boost.transform(D0.momentum()).p3().unit();
	double ctheta = e1z.dot(axis1);
	_h_ctheta.fill(xE,ctheta);
	double phi = atan2(e1y.dot(axis1),e1x.dot(axis1));
	_h_phi.fill(xE,phi);
	_h_01.fill(xE,ctheta,cos(phi));
	
	if(xE>0.3) {
	  _h_ctheta_all->fill(ctheta);
	  _h_phi_all->fill(phi);
	  _h_01_all->fill(ctheta,cos(phi));
	}
      }
    }

    pair<double,double> calcRho(Histo1DPtr hist,unsigned int mode) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai,bi;
	if(mode==0) {
	  ai = 0.25*(bin.xMax()*(3.-sqr(bin.xMax())) - bin.xMin()*(3.-sqr(bin.xMin())));
	  bi = 0.75*(bin.xMin()*(1.-sqr(bin.xMin())) - bin.xMax()*(1.-sqr(bin.xMax())));
	}
	else if(mode==1) {
	  ai = 0.5/M_PI*(bin.xMax()-bin.xMin());
	  bi = 0.5/M_PI*(sin(2.*bin.xMin())-sin(2.*bin.xMax()));
	}
	else  {
	  ai=0.;
	  bi= sqrt(0.5)*((1.-sqr(bin.xMax()))*sqrt(1.-sqr(bin.xMax()))-
			 (1.-sqr(bin.xMin()))*sqrt(1.-sqr(bin.xMin())));
	}
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    pair<double,pair<double,double> > calcAlpha(Histo1DPtr hist) {
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
    
    /// Normalise histograms etc., after the run
    void finalize() {
      vector<double> x = {0.3,0.4,0.5,0.6,0.7,0.8,1.};
      Scatter2DPtr h_alpha, h_rho00, h_rhooff, h_01;
      book(h_alpha , 1,1,1);
      book(h_rho00 , 2,1,1);
      book(h_rhooff, 2,1,2);
      book(h_01    , 2,1,3);
      for(unsigned int ix=0;ix<6;++ix) {
	double integral = _h_ctheta.histos()[ix]->integral();
	scale( _h_ctheta.histos()[ix],1./integral);
	scale( _h_01.histos()[ix],1./integral);
	normalize(_h_phi.histos()[ix]);
	pair<double,double> rho00 = calcRho(_h_ctheta.histos()[ix],0);
	h_rho00->addPoint(0.5*(x[ix]+x[ix+1]), rho00.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			  make_pair(rho00.second,rho00.second) );
	pair<double,pair<double,double> > alpha = calcAlpha(_h_ctheta.histos()[ix]);
	h_alpha->addPoint(0.5*(x[ix]+x[ix+1]), alpha.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			  alpha.second);
	pair<double,double> rhooff = calcRho(_h_phi.histos()[ix],1);
	h_rhooff->addPoint(0.5*(x[ix]+x[ix+1]), rhooff.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			  make_pair(rhooff.second,rhooff.second) );
	pair<double,double> rho01 = calcRho(_h_01.histos()[ix],2);
	h_01->addPoint(0.5*(x[ix]+x[ix+1]), rho01.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
		       make_pair(rho01.second,rho01.second) );
      }
      // integral over z
      book(h_alpha , 1,2,1);
      book(h_rho00 , 2,2,1);
      book(h_rhooff, 2,2,2);
      book(h_01    , 2,2,3);
      double integral = _h_ctheta_all->integral();
      scale( _h_ctheta_all,1./integral);
      scale( _h_01_all,1./integral);
      normalize(_h_phi_all);
      pair<double,double> rho00 = calcRho(_h_ctheta_all,0);
      h_rho00->addPoint(.65, rho00.first, make_pair(0.35,0.35),
			make_pair(rho00.second,rho00.second) );
      pair<double,pair<double,double> > alpha = calcAlpha(_h_ctheta_all);
      h_alpha->addPoint(.65, alpha.first, make_pair(0.35,0.35),
			alpha.second);
      pair<double,double> rhooff = calcRho(_h_phi_all,1);
      h_rhooff->addPoint(.65, rhooff.first, make_pair(0.35,0.35),
			 make_pair(rhooff.second,rhooff.second) );
      pair<double,double> rho01 = calcRho(_h_01_all,2);
      h_01->addPoint(.65, rho01.first, make_pair(0.35,0.35),
		     make_pair(rho01.second,rho01.second) );
    }

    //@}


    /// @name Histograms
    //@{
    BinnedHistogram _h_ctheta,_h_phi,_h_01;
    Histo1DPtr _h_ctheta_all,_h_phi_all,_h_01_all;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TPC_1991_I316132);


}
