// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D* polarization in B0 -> D* tau nu_tau
  class BELLE_2019_I1724068 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2019_I1724068);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_cTheta,1,1,1);
    }

    void findChildren(const Particle & p, unsigned int & nprod,
		      Particles & Dstar, Particles & tau, Particles & nu) {
      for(const Particle & child : p.children()) {
	if(child.pid()==-413) {
	  ++nprod;
	  Dstar.push_back(child);
	}
	else if(child.pid()==-15) {
	  ++nprod;
	  tau.push_back(child);
	}
	else if(child.pid()==16) {
	  ++nprod;
	  nu.push_back(child);
	}
	else if(child.pid()==22)
	  continue;
	else if(child.children().empty() ||
		child.pid()==111 || child.pid()==221 || child.pid()==331) {
	  ++nprod;
	}
	else {
	  findChildren(child,nprod,Dstar,tau,nu);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over B0 mesons 
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::B0)) {
	// find the B decay
	unsigned int nprod = 0;
	Particles Dstar,tau,nu;
	findChildren(p,nprod,Dstar,tau,nu);
	if(nprod!=3 || Dstar.size()!=1 || tau.size() !=1 || nu.size()!=1)
	  continue;
	// and the D* decay
	if(Dstar[0].children().size()!=2) continue;
	Particle D0;
	if(Dstar[0].children()[0].pid()==-421 &&
	   Dstar[0].children()[1].pid()==-211) {
	  D0 = Dstar[0].children()[0];
	}
	else if(Dstar[0].children()[1].pid()==-421 &&
		Dstar[0].children()[0].pid()==-211) {
	  D0 = Dstar[0].children()[1];
	}
	else
	  continue;
	// compute the helicity angle
	// boost to B0 rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum pDstar = boost1.transform(Dstar[0].momentum());
	FourMomentum pD     = boost1.transform(D0      .momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pDstar.betaVec());
	pD = boost2.transform(pD);
	double cTheta = pD.p3().unit().dot(pDstar.p3().unit());
	if(cTheta<=0.) _h_cTheta->fill(cTheta);
      }
    }

    pair<double,double> calcF(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
      	double Oi = bin.area();
      	if(Oi==0.) continue;
      	double ai = 0.5*(bin.xMin()*(sqr(bin.xMin())-3.)-bin.xMax()*(sqr(bin.xMax())-3.));
	double bi = 1.5*(bin.xMin()*(1.-sqr(bin.xMin()))-
			 bin.xMax()*(1.-sqr(bin.xMax())));
      	double Ei = bin.areaErr();
       	sum1 += sqr(bi/Ei);
      	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_cTheta);
      Scatter2DPtr _h_F;
      book(_h_F,2,1,1);
      pair<double,double> F = calcF(_h_cTheta);
      _h_F->addPoint(0.5, F.first, make_pair(0.5,0.5), make_pair(F.second,F.second) );
      
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_cTheta;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2019_I1724068);

}
