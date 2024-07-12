// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Xi_c0 -> Xi-pi+ asymmetry
  class CLEO_2000_I537236 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2000_I537236);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      // Book histograms
      book(_h_ctheta, "ctheta"  , 20,-1,1);
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
	_h_ctheta->fill(cTheta,1.);
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
      normalize(_h_ctheta);
      Scatter2DPtr _h_alpha;
      book(_h_alpha,1,1,1);
      pair<double,double> alpha = calcAlpha(_h_ctheta);
      _h_alpha->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_ctheta;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(CLEO_2000_I537236);


}
