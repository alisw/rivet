// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 decay asymmetries
  class BELLE_2021_I1859517 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1859517);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // book histograms
      book(_h_Lambda,1,1,1);
      book(_h_Sigma0,2,1,1);
      book(_h_Sigmap,3,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Xi_c0 baryons
      for( const Particle& Xic : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4132)) {
	int sign = Xic.pid()/4132;
	if(Xic.children().size()!=2) continue;
	Particle baryon1,meson1;
	if(Xic.children()[0].pid()==sign*3122 && 
	   Xic.children()[1].pid()==-sign*313) {
	  baryon1 = Xic.children()[0];
	  meson1  = Xic.children()[1];
	}
	else if(Xic.children()[1].pid()==sign*3122 && 
		Xic.children()[0].pid()==-sign*313) {
	  baryon1 = Xic.children()[1];
	  meson1  = Xic.children()[0];
	}
	else if(Xic.children()[0].pid()==sign*3212 && 
		Xic.children()[1].pid()==-sign*313) {
	  baryon1 = Xic.children()[0];
	  meson1  = Xic.children()[1];
	}
	else if(Xic.children()[1].pid()==sign*3212 && 
		Xic.children()[0].pid()==-sign*313) {
	  baryon1 = Xic.children()[1];
	  meson1  = Xic.children()[0];
	}
	else if(Xic.children()[0].pid()==sign*3222 && 
		Xic.children()[1].pid()==-sign*323) {
	  baryon1 = Xic.children()[0];
	  meson1  = Xic.children()[1];
	}
	else if(Xic.children()[1].pid()==sign*3222 && 
		Xic.children()[0].pid()==-sign*323) {
	  baryon1 = Xic.children()[1];
	  meson1  = Xic.children()[0];
	}
	else
	  continue;
	Particle baryon2,meson2;
	if(baryon1.abspid()==3122) {
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
	else if(baryon1.abspid()==3212) {
	  if(baryon1.children()[0].pid()== sign*3122 && 
	     baryon1.children()[1].pid()==22) {
	    baryon2 = baryon1.children()[0];
	    meson2  = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()== sign*3122 && 
		  baryon1.children()[0].pid()==22) {
	    baryon2 = baryon1.children()[1];
	    meson2  = baryon1.children()[0];
	  }
	  else
	    continue;
	}
	else if(baryon1.abspid()==3222) {
	  if(baryon1.children()[0].pid()== sign*2212 && 
	     baryon1.children()[1].pid()==111) {
	    baryon2 = baryon1.children()[0];
	    meson2  = baryon1.children()[1];
	  }
	  else if(baryon1.children()[1].pid()== sign*2212 && 
		  baryon1.children()[0].pid()==111) {
	    baryon2 = baryon1.children()[1];
	    meson2  = baryon1.children()[0];
	  }
	  else
	    continue;
	}
	// first boost to the Xic rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Xic.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(baryon1.abspid()==3122) {
	  _h_Lambda->fill(cTheta);
	}
	else if(baryon1.abspid()==3212) {
	  _h_Sigma0->fill(cTheta);
	}
	else if(baryon1.abspid()==3222) {
	  _h_Sigmap->fill(cTheta);
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

    /// Normalise histograms etc., after the run
    void finalize() {
      // Lambda Kbar*0
      normalize(_h_Lambda);
      Scatter2DPtr _h_alpha;
      book(_h_alpha,4,1,1);
      pair<double,double> alpha = calcAlpha(_h_Lambda);
      _h_alpha->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      // Sigma0 Kbar*0
      normalize(_h_Sigma0);
      book(_h_alpha,4,1,2);
      alpha = calcAlpha(_h_Sigma0);
      _h_alpha->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      // Sigma+ K*-
      normalize(_h_Sigmap);
      book(_h_alpha,4,1,3);
      alpha = calcAlpha(_h_Sigmap);
      _h_alpha->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );

    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Lambda,_h_Sigma0,_h_Sigmap;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1859517);

}
