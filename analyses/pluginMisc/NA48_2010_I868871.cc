// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include <sstream>
namespace Rivet {


  /// @brief asymmetrics in Xi0 decays
  class NA48_2010_I868871 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(NA48_2010_I868871);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // Book histograms
      book(_h_ctheta_pi0  , "ctheta_pi0"  , 20,-1,1);
      book(_h_ctheta_gamma, "ctheta_gamma", 20,-1,1);
      double step=0.1;
      double xmin=-1.;
      for(unsigned int ix=0;ix<20;++ix) {
	Histo1DPtr temp;
	std::ostringstream title;
	title << "ctheta_Sigma_" << ix;
	book(temp,title.str(), 20,-1,1);
	_h_ctheta_Sigma.add(xmin, xmin+step, temp);
	xmin+=step;
      }
      _nSigma=0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Omega baryons
      for(const Particle& Xi : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==3322)) {
	int sign = Xi.pid()/3322;
	if(Xi.children().size()!=2) continue;
	Particle baryon1,meson1;
	unsigned int mode(0);
	if(Xi.children()[0].pid()==sign*3122 && 
	   Xi.children()[1].pid()==111) {
	  baryon1 = Xi.children()[0];
	  meson1  = Xi.children()[1];
	  mode=1;
	}
	else if(Xi.children()[1].pid()==sign*3122 && 
		Xi.children()[0].pid()==111) {
	  baryon1 = Xi.children()[1];
	  meson1  = Xi.children()[0];
	  mode=1;
	}
	else if(Xi.children()[0].pid()==sign*3122 && 
	   Xi.children()[1].pid()==22) {
	  baryon1 = Xi.children()[0];
	  meson1  = Xi.children()[1];
	  mode=2;
	}
	else if(Xi.children()[1].pid()==sign*3122 && 
		Xi.children()[0].pid()==22) {
	  baryon1 = Xi.children()[1];
	  meson1  = Xi.children()[0];
	  mode=2;
	}
	else if(Xi.children()[0].pid()==sign*3212 && 
	   Xi.children()[1].pid()==22) {
	  baryon1 = Xi.children()[0];
	  meson1  = Xi.children()[1];
	  mode=3;
	}
	else if(Xi.children()[1].pid()==sign*3212 && 
		Xi.children()[0].pid()==22) {
	  baryon1 = Xi.children()[1];
	  meson1  = Xi.children()[0];
	  mode=3;
	}
	else
	  continue;
	if(baryon1.children().size()!=2) continue;
	Particle baryon2,meson2,baryon3,meson3;
	if(mode==1 || mode ==2) {
	  if(baryon1.children()[0].pid()== sign*2212 && 
	     baryon1.children()[1].pid()==-sign*211) {
	    baryon2 = baryon1.children()[0];
	    meson2  = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()== sign*2212 && 
		  baryon1.children()[0].pid()==-sign*211) {
	    baryon2 = baryon1.children()[1];
	    meson2  = baryon1.children()[0];
	  }
	  else
	    continue;
	}
	else if(mode==3) {
	  if(baryon1.children()[0].pid()== sign*3122 && 
	     baryon1.children()[1].pid()== 22) {
	    baryon2 = baryon1.children()[0];
	    meson2  = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()== sign*3122 && 
		  baryon1.children()[0].pid()== 22) {
	    baryon2 = baryon1.children()[1];
	    meson2  = baryon1.children()[0];
	  }
	  else
	    continue;
	  if(baryon2.children()[0].pid()== sign*2212 && 
	     baryon2.children()[1].pid()==-sign*211) {
	    baryon3 = baryon2.children()[0];
	    meson3  = baryon2.children()[1];
	  }
	  else if(baryon2.children()[1].pid()== sign*2212 && 
		  baryon2.children()[0].pid()==-sign*211) {
	    baryon3 = baryon2.children()[1];
	    meson3  = baryon2.children()[0];
	  }
	  else
	    continue;
	}
	// first boost to the Xi rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Xi.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(mode==1) {
	  _h_ctheta_pi0->fill(cTheta,1.);
	}
	else if(mode==2) {
	  _h_ctheta_gamma->fill(cTheta,1.);
	}
	else if(mode==3) {
	  FourMomentum pbaryon3 = boost1.transform(baryon3.momentum());
	  FourMomentum pp2      = boost2.transform(pbaryon3);
	  Vector3 axis2 = pp.p3().unit();
	  double cTheta2 = pp2.p3().unit().dot(axis2);
	  _h_ctheta_Sigma.fill(cTheta,cTheta2,1.);
	  _nSigma += 1.;
	}
      }
    }

    pair<double,double> calcAlpha(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai = 0.5*(bin.xMax()-bin.xMin());
	double bi = 0.5*ai*(bin.xMax()+bin.xMin());
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }
    
    pair<double,double> calcAlpha(BinnedHistogram hist) {
      double sum1(0.),sum2(0.);
      double step=0.1;
      double xmin=-1.;
      for(unsigned int ix=0;ix<20;++ix) {
	double xsum=2.*xmin+step;
	Histo1DPtr h2 = hist.histo(xmin+0.5*step);
	for (auto bin : h2->bins() ) {
	  double Oi = bin.area();
	  if(Oi==0.) continue;
	  double ai = 0.25*(bin.xMax()-bin.xMin())*step;
	  double bi = 0.25*ai*(bin.xMax()+bin.xMin())*xsum;
	  double Ei = bin.areaErr();
	  sum1 += sqr(bi/Ei);
	  sum2 += bi/sqr(Ei)*(Oi-ai);
	}
	xmin+=step;
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Xi0 -> Lambda0 pi0
      normalize(_h_ctheta_pi0);
      Scatter2DPtr _h_alpha_pi0;
      book(_h_alpha_pi0,1,1,1);
      pair<double,double> alpha = calcAlpha(_h_ctheta_pi0);
      _h_alpha_pi0->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      // Xi0 -> Lambda gamma (N.B. sign due defns)
      normalize(_h_ctheta_gamma);
      Scatter2DPtr _h_alpha_gamma;
      book(_h_alpha_gamma,2,1,1);
      alpha = calcAlpha(_h_ctheta_gamma);
      _h_alpha_gamma->addPoint(0.5,-alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      // Xi0 -> Sigma gamma
      _h_ctheta_Sigma.scale(1./_nSigma,this);
      Scatter2DPtr _h_alpha_Sigma;
      book(_h_alpha_Sigma,3,1,1);
      alpha = calcAlpha(_h_ctheta_Sigma);
      _h_alpha_Sigma->addPoint(0.5,alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ctheta_pi0,_h_ctheta_gamma;
    BinnedHistogram _h_ctheta_Sigma;
    double _nSigma;
    //@}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(NA48_2010_I868871);


}
