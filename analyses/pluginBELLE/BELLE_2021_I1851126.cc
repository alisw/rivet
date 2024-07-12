// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 -> Xi-pi+ asymmetry
  class BELLE_2021_I1851126 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2021_I1851126);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // Book histograms
      book(_h_c_P,1,1,1);
      book(_h_c_M,1,1,2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Xi_c0 baryons
      for( const Particle& Xic : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4132)) {
	int sign = Xic.pid()/4132;
	if(Xic.children().size()!=2) continue;
	Particle baryon1,meson1;
	if(Xic.children()[0].pid()==sign*3312 && 
	   Xic.children()[1].pid()==sign*211) {
	  baryon1 = Xic.children()[0];
	  meson1  = Xic.children()[1];
	}
	else if(Xic.children()[1].pid()==sign*3312 && 
		Xic.children()[0].pid()==sign*211) {
	  baryon1 = Xic.children()[1];
	  meson1  = Xic.children()[0];
	}
	else
	  continue;
	Particle baryon2,meson2;
	if(baryon1.children()[0].pid()== sign*3122 && 
	   baryon1.children()[1].pid()==-sign*211) {
	  baryon2 = baryon1.children()[0];
	  meson2  = baryon1.children()[1];
	}
	else if(baryon1.children()[1].pid()== sign*3122 && 
		baryon1.children()[0].pid()==-sign*211) {
	  baryon2 = baryon1.children()[1];
	  meson2  = baryon1.children()[0];
	}
	else
	  continue;
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
	if(sign>0)
	  _h_c_P->fill(cTheta,1.);
	else
	  _h_c_M->fill(cTheta,1.);
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
      // first mode
      normalize(_h_c_P);
      Scatter2DPtr _h_alpha_P;
      book(_h_alpha_P,2,1,1);
      pair<double,double> alphaP = calcAlpha(_h_c_P);
      alphaP.first /= -0.401;
      alphaP.second/= -0.401;
      _h_alpha_P->addPoint(0.5, alphaP.first, make_pair(0.5,0.5), make_pair(alphaP.second,alphaP.second) );
      // second mode
      normalize(_h_c_M);
      Scatter2DPtr _h_alpha_M;
      book(_h_alpha_M,2,1,2);
      pair<double,double> alphaM = calcAlpha(_h_c_M);
      alphaM.first /= 0.389;
      alphaM.second/= 0.389;
      _h_alpha_M->addPoint(0.5, alphaM.first, make_pair(0.5,0.5), make_pair(alphaM.second,alphaM.second) );
      // average
      double aver = 0.5*(-alphaP.first+alphaM.first);
      double err  = 0.5*sqrt(sqr(alphaP.second)+sqr(alphaM.second));
      Scatter2DPtr _h_alpha_aver;
      book(_h_alpha_aver,2,1,3);
      _h_alpha_aver->addPoint(0.5, aver, make_pair(0.5,0.5), make_pair(err,err) );
      // asymetry
      double asym = (alphaP.first+alphaM.first)/(alphaP.first-alphaM.first);
      err         = 2./sqr(alphaP.first-alphaM.first)*sqrt(sqr(alphaM.first *alphaP.second)+
							   sqr(alphaM.second*alphaP.first ));
      Scatter2DPtr _h_alpha_asym;
      book(_h_alpha_asym,2,1,4);
      _h_alpha_asym->addPoint(0.5, asym, make_pair(0.5,0.5), make_pair(err,err) );
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_c_M,_h_c_P;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2021_I1851126);

}
