// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda_c -> Sigma+ pi0,eta,eta' decay asymmetries
  class BELLE_2022_I2140379 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2140379);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // histograms
      for(unsigned int ix=0;ix<3;++ix)
	book(_h[ix],2,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_c baryons
      for( const Particle& Lambdac : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4122)) {
	int sign = Lambdac.pid()/4122;
	if(Lambdac.children().size()!=2) continue;
	Particle baryon1;
	int imeson=-1;
	if(Lambdac.children()[0].pid()==sign*3222 && 
	   Lambdac.children()[1].pid()==111) {
	  baryon1 = Lambdac.children()[0];
	  imeson=0;
	}
	else if(Lambdac.children()[1].pid()==sign*3222 && 
		Lambdac.children()[0].pid()==111) {
	  baryon1 = Lambdac.children()[1];
	  imeson=0;
	}
	else if(Lambdac.children()[0].pid()==sign*3222 && 
		Lambdac.children()[1].pid()==221) {
	  baryon1 = Lambdac.children()[0];
	  imeson=1;
	}
	else if(Lambdac.children()[1].pid()==sign*3222 && 
		Lambdac.children()[0].pid()==221) {
	  baryon1 = Lambdac.children()[1];
	  imeson=1;
	}
	else if(Lambdac.children()[0].pid()==sign*3222 && 
		Lambdac.children()[1].pid()==331) {
	  baryon1 = Lambdac.children()[0];
	  imeson=2;
	}
	else if(Lambdac.children()[1].pid()==sign*3222 && 
		Lambdac.children()[0].pid()==331) {
	  baryon1 = Lambdac.children()[1];
	  imeson=2;
	}
	else
	  continue;
	Particle baryon2;
	if(baryon1.children()[0].pid()== sign*2212 && 
	   baryon1.children()[1].pid()== 111) {
	  baryon2 = baryon1.children()[0];
	}
	else if(baryon1.children()[1].pid()== sign*2212 && 
		baryon1.children()[0].pid()== 111) {
	  baryon2 = baryon1.children()[1];
	}
	else
	  continue;
	// first boost to the Lambdac rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambdac.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(baryon1.momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2.momentum());
	// to sigma+ rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h[imeson]->fill(cTheta);
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
      pair<double,double> aSigma(-.983,0.013);
      for(unsigned int ix=0;ix<3;++ix) {
	normalize(_h[ix]);
	Scatter2DPtr _h_alpha1;
	book(_h_alpha1,1,1+ix,1);
	pair<double,double> alpha = calcAlpha(_h[ix]);
	_h_alpha1->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
	// divide out alpha Sigma
	alpha.second = alpha.first/aSigma.first*
	  sqrt(sqr(alpha.second/alpha.first) + sqr(aSigma.second/aSigma.first));
	alpha.first /= aSigma.first;
	Scatter2DPtr _h_alpha2;
	book(_h_alpha2,1,1+ix,2);
	_h_alpha2->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2140379);

}
