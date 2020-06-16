// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi-> Lambda pi asymmetry
  class E756_2000_I530367 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(E756_2000_I530367);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );

      // Book histograms
      book(_h_cthetaP, "cthetaP",20,-1,1);
      book(_h_cthetaM, "cthetaM",20,-1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Xi- baryons
      for (const Particle& Xi : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==3312)) {
	int sign = Xi.pid()/3312;
	if(Xi.children().size()!=2) continue;
	Particle Lambda,pion1;
	if(Xi.children()[0].pid()==sign*3122 && 
	   Xi.children()[1].pid()==-sign*211) {
	  Lambda = Xi.children()[0];
	  pion1   = Xi.children()[1];
	}
	else if(Xi.children()[1].pid()==sign*3122 && 
		Xi.children()[0].pid()==-sign*211) {
	  Lambda = Xi.children()[1];
	  pion1   = Xi.children()[0];
	}
	else
	  continue;
	if(Lambda.children().size()!=2) continue;
	Particle proton,pion2;
	if(Lambda.children()[0].pid()==sign*2212 && 
	   Lambda.children()[1].pid()==-sign*211) {
	  proton = Lambda.children()[0];
	  pion2   = Lambda.children()[1];
	}
	else if(Lambda.children()[1].pid()==sign*2212 && 
		Lambda.children()[0].pid()==-sign*211) {
	  proton = Lambda.children()[1];
	  pion2   = Lambda.children()[0];
	}
	else
	  continue;
	// boost to xi rest frame first
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Xi.momentum().betaVec());
	FourMomentum pLambda = boost1.transform(Lambda.momentum());
	FourMomentum pproton = boost1.transform(proton.momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pLambda.betaVec());
	Vector3 axis = pLambda.p3().unit();
	FourMomentum pp = boost2.transform(pproton);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	if(sign==1) {
	  _h_cthetaM->fill(cTheta);
	}
	else {
	  _h_cthetaP->fill(cTheta);
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
      normalize(_h_cthetaP);
      normalize(_h_cthetaM);
      // calculate the values of alpha
      // xibar+
      Scatter2DPtr _h_alphaP;
      book(_h_alphaP, 1,1,2);
      pair<double,double> alpha = calcAlpha(_h_cthetaP);
      _h_alphaP->addPoint(0.5, alpha.first, make_pair(0.5,0.5),
			  make_pair(alpha.second,alpha.second) );
      // xi-
      Scatter2DPtr _h_alphaM;
      book(_h_alphaM, 1,1,1);
      alpha = calcAlpha(_h_cthetaM);
      _h_alphaM->addPoint(0.5, alpha.first, make_pair(0.5,0.5),
			  make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_cthetaP,_h_cthetaM;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(E756_2000_I530367);


}
