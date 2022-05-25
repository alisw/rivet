// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief asymmetry in Omega-> Lambda K
  class HyperCP_2005_I677384 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HyperCP_2005_I677384);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );

      // Book histograms
      book(_h_cthetaP  , "cthetaP"  ,20,-1,1);
      book(_h_cthetaM  , "cthetaM"  ,20,-1,1);
      book(_h_cthetaAll, "cthetaAll",20,-1,1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Omega baryons
      for(const Particle& Omega : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==3334)) {
	int sign = Omega.pid()/3334;
	if(Omega.children().size()!=2) continue;
	Particle Lambda,kaon;
	if(Omega.children()[0].pid()==sign*3122 && 
	   Omega.children()[1].pid()==-sign*321) {
	  Lambda = Omega.children()[0];
	  kaon   = Omega.children()[1];
	}
	else if(Omega.children()[1].pid()==sign*3122 && 
		Omega.children()[0].pid()==-sign*321) {
	  Lambda = Omega.children()[1];
	  kaon   = Omega.children()[0];
	}
	else
	  continue;
	if(Lambda.children().size()!=2) continue;
	Particle proton,pion;
	if(Lambda.children()[0].pid()==sign*2212 && 
	   Lambda.children()[1].pid()==-sign*211) {
	  proton = Lambda.children()[0];
	  pion   = Lambda.children()[1];
	}
	else if(Lambda.children()[1].pid()==sign*2212 && 
		Lambda.children()[0].pid()==-sign*211) {
	  proton = Lambda.children()[1];
	  pion   = Lambda.children()[0];
	}
	else
	  continue;
	// first boost to the Omega rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Omega.momentum().betaVec());
	FourMomentum pLambda = boost1.transform(Lambda.momentum());
	FourMomentum pproton = boost1.transform(proton.momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pLambda.betaVec());
	Vector3 axis = pLambda.p3().unit();
	FourMomentum pp = boost2.transform(pproton);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h_cthetaAll->fill(cTheta,1.);
	if(sign==1) {
	  _h_cthetaM->fill(cTheta,1.);
	}
	else {
	  _h_cthetaP->fill(cTheta,1.);
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
      normalize(_h_cthetaP  );
      normalize(_h_cthetaM  );
      normalize(_h_cthetaAll);
      // calculate the values of alpha
      Scatter2DPtr _h_alphaP;
      book(_h_alphaP,1,1,1);
      pair<double,double> alpha = calcAlpha(_h_cthetaP);
      _h_alphaP->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      Scatter2DPtr _h_alphaM;
      book(_h_alphaM,1,1,2);
      alpha = calcAlpha(_h_cthetaM);
      _h_alphaM->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
      Scatter2DPtr _h_alphaAll;
      book(_h_alphaAll,1,1,3);
      alpha = calcAlpha(_h_cthetaAll);
      _h_alphaAll->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_cthetaP,_h_cthetaM,_h_cthetaAll;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(HyperCP_2005_I677384);


}
