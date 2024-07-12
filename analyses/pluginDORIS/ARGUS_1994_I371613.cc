// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  Lambda_c -> Lambda l+ nu_l asymmetry
  class ARGUS_1994_I371613 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1994_I371613);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS" );
      
      // Book histograms
      book(_h_Lambda, "/TMP/hLambda", 20,-1.,1.);

    }

    void findChildren(Particle parent, int sign, unsigned int & npart,
		      Particles & lambda, Particles & e, Particles & nu) {
      for(const Particle & child : parent.children()) {
	if(child.pid()==sign*PID::LAMBDA) {
	  lambda.push_back(child);
	  ++npart;
	}
	else if(child.pid()==-sign*PID::EMINUS || child.pid()==-sign*PID::MUON) {
	  e.push_back(child);
	  ++npart;
	}
	else if(child.pid()==sign*PID::NU_E || child.pid()==sign*PID::NU_MU) {
	  nu.push_back(child);
	  ++npart;
	}
	else if(!child.children().empty()) {
	  findChildren(child,sign,npart,lambda,e,nu);
	}
	else {
	  ++npart;
	}
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // loop over Lambda_c baryons
      for( const Particle& Lambdac : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==4122)) {
	int sign = Lambdac.pid()/4122;
        Particles lambda,e,nu;
        unsigned int npart(0);
        findChildren(Lambdac,sign,npart,lambda,e,nu);
	if(npart!=3 || lambda.size()!=1 || e.size()!=1 || nu.size()!=1) continue;
	Particle baryon2;
	if(lambda[0].children()[0].pid()== sign*2212 && 
	   lambda[0].children()[1].pid()== -sign*211) {
	  baryon2 = lambda[0].children()[0];
	}
	else if(lambda[0].children()[1].pid()== sign*2212 && 
		lambda[0].children()[0].pid()== -sign*211) {
	  baryon2 = lambda[0].children()[1];
	}
	else
	  continue;
	// mass cut
	double mLL = (lambda[0].momentum()+e[0].momentum()).mass();
	if(mLL<1.85 || mLL>2.2) continue;
	// first boost to the Lambdac rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambdac.momentum().betaVec());
	FourMomentum pbaryon1 = boost1.transform(lambda[0].momentum());
	FourMomentum pbaryon2 = boost1.transform(baryon2  .momentum());
	// to lambda rest frame
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pbaryon1.betaVec());
	Vector3 axis = pbaryon1.p3().unit();
	FourMomentum pp = boost2.transform(pbaryon2);
	// calculate angle
	double cTheta = pp.p3().unit().dot(axis);
	_h_Lambda->fill(cTheta);
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
      //  asymmetry
      normalize(_h_Lambda);
      Scatter2DPtr _h_alpha;
      book(_h_alpha,1,1,1);
      pair<double,double> alpha = calcAlpha(_h_Lambda);
      _h_alpha->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_Lambda;
    //@}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1994_I371613);

}
