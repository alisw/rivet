// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c -> Lambda pi and Lambda_c _> Sigma+ pi0 asymmetries
  class CLEO_1995_I392704 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1995_I392704);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      
      // Book histograms
      book(_h_Lambda, 1,1,1);
      book(_h_Sigma , 2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_c baryons
      for( const Particle& Lambdac : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4122)) {
	int sign = Lambdac.pid()/4122;
	if(Lambdac.children().size()!=2) continue;
	Particle baryon1;
	bool lambda=true;
	if(Lambdac.children()[0].pid()==sign*3122 && 
	   Lambdac.children()[1].pid()==sign*211) {
	  baryon1 = Lambdac.children()[0];
	}
	else if(Lambdac.children()[1].pid()==sign*3122 && 
		Lambdac.children()[0].pid()==sign*211) {
	  baryon1 = Lambdac.children()[1];
	}
	else if(Lambdac.children()[0].pid()==sign*3222 && 
		Lambdac.children()[1].pid()==111) {
	  baryon1 = Lambdac.children()[0];
	  lambda=false;
	}
	else if(Lambdac.children()[1].pid()==sign*3222 && 
		Lambdac.children()[0].pid()==111) {
	  baryon1 = Lambdac.children()[0];
	  lambda=false;
	}
	else
	  continue;
	int idMeson = lambda ? -sign*211 : 111;
	Particle baryon2;
	if(baryon1.children()[0].pid()== sign*2212 && 
	   baryon1.children()[1].pid()== idMeson) {
	  baryon2 = baryon1.children()[0];
	}
	else if(baryon1.children()[1].pid()== sign*2212 && 
		baryon1.children()[0].pid()== idMeson) {
	  baryon2 = baryon1.children()[1];
	}
	else
	  continue;
	// first boost to the Lambdac rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambdac.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(lambda)
	  _h_Lambda->fill(cTheta);
	else
	  _h_Sigma->fill(cTheta);
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
      // Lambda_c -> Lambda pi+
      normalize(_h_Lambda);
      Scatter2DPtr _h_alpha1;
      book(_h_alpha1,3,1,1);
      pair<double,double> alpha = calcAlpha(_h_Lambda);
      _h_alpha1->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      // Lambda_c -> Sigma+ pi0
      normalize(_h_Sigma);
      Scatter2DPtr _h_alpha2;
      book(_h_alpha2,4,1,1);
      alpha = calcAlpha(_h_Sigma);
      _h_alpha2->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Lambda, _h_Sigma;
    //@}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1995_I392704);

}
