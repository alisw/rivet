// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Omega decay asymmetries
  class WA46_1984_I206647 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(WA46_1984_I206647);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );

      // Book histograms
      book(_h_cthetalam, "cthetaLambda",20,-1,1);
      book(_h_cthetaxi0, "cthetaXi0"   ,20,-1,1);
      book(_h_cthetaxim, "cthetaXim"   ,20,-1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Omega baryons
      for(const Particle& Omega : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==3334)) {
	int sign = Omega.pid()/3334;
	if(Omega.children().size()!=2) continue;
	Particle baryon1,meson1;
	if(Omega.children()[0].pid()==sign*3122 && 
	   Omega.children()[1].pid()==-sign*321) {
	  baryon1 = Omega.children()[0];
	  meson1 = Omega.children()[1];
	}
	else if(Omega.children()[1].pid()==sign*3122 && 
		Omega.children()[0].pid()==-sign*321) {
	  baryon1 = Omega.children()[1];
	  meson1 = Omega.children()[0];
	}
	else if(Omega.children()[0].pid()==sign*3322 && 
		Omega.children()[1].pid()==-sign*211) {
	  baryon1 = Omega.children()[0];
	  meson1 = Omega.children()[1];
	}
	else if(Omega.children()[1].pid()==sign*3322 && 
		Omega.children()[0].pid()==-sign*211) {
	  baryon1 = Omega.children()[1];
	  meson1 = Omega.children()[0];
	}
	else if(Omega.children()[0].pid()==sign*3312 && 
		Omega.children()[1].pid()==111) {
	  baryon1 = Omega.children()[0];
	  meson1 = Omega.children()[1];
	}
	else if(Omega.children()[1].pid()==sign*3312 && 
		Omega.children()[0].pid()==111) {
	  baryon1 = Omega.children()[1];
	  meson1 = Omega.children()[0];
	}
	else
	  continue;
	if(baryon1.children().size()!=2) continue;
	Particle baryon2,meson2;
	if(baryon1.abspid()==3122) {
	  if(baryon1.children()[0].pid()==sign*2212 && 
	     baryon1.children()[1].pid()==-sign*211) {
	    baryon2 = baryon1.children()[0];
	    meson2   = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()==sign*2212 && 
		  baryon1.children()[0].pid()==-sign*211) {
	    baryon2 = baryon1.children()[1];
	    meson2   = baryon1.children()[0];
	  }
	  else
	    continue;
	}
	else if(baryon1.abspid()==3322) {
	  if(baryon1.children()[0].pid()==sign*3122 && 
	     baryon1.children()[1].pid()==111) {
	    baryon2 = baryon1.children()[0];
	    meson2   = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()==sign*3122 && 
		  baryon1.children()[0].pid()==111) {
	    baryon2 = baryon1.children()[1];
	    meson2   = baryon1.children()[0];
	  }
	  else
	    continue;
	}
	else if (baryon1.abspid()==3312) {
	  if(baryon1.children()[0].pid()==sign*3122 && 
	     baryon1.children()[1].pid()==-sign*211) {
	    baryon2 = baryon1.children()[0];
	    meson2   = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()==sign*3122 && 
		  baryon1.children()[0].pid()==-sign*211) {
	    baryon2 = baryon1.children()[1];
	    meson2   = baryon1.children()[0];
	  }
	  else
	    continue;
	}
	// first boost to the Omega rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Omega.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(baryon1.abspid()==3122)
	  _h_cthetalam->fill(cTheta);
	else if(baryon1.abspid()==3322)
	  _h_cthetaxi0->fill(cTheta);
	else if(baryon1.abspid()==3312)
	  _h_cthetaxim->fill(cTheta);
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
    
    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_cthetalam);
      normalize(_h_cthetaxi0);
      normalize(_h_cthetaxim);
      // calculate the values of alpha
      Scatter2DPtr _h_alphaLam;
      book(_h_alphaLam,1,1,1);
      pair<double,double> alpha = calcAlpha(_h_cthetalam);
      _h_alphaLam->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      Scatter2DPtr _h_alphaXi0;
      book(_h_alphaXi0,1,1,2);
      alpha = calcAlpha(_h_cthetaxi0);
      _h_alphaXi0->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      Scatter2DPtr _h_alphaXim;
      book(_h_alphaXim,1,1,3);
      alpha = calcAlpha(_h_cthetaxim);
      _h_alphaXim->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_cthetalam,_h_cthetaxi0,_h_cthetaxim;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(WA46_1984_I206647);


}
