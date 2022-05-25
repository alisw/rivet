// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau polarization in B -> D* tau nu_tau
  class BELLE_2018_I1621272 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2018_I1621272);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_pi ,"/TMP/PI" ,20,-1.,1.);
      book(_h_rho,"/TMP/RHO",20,-1.,1.);
    }

    void findChildren(const Particle & p, int & sign, unsigned int & nprod,
		      Particles & Dstar, Particles & tau, Particles & nu) {
      for(const Particle & child : p.children()) {
	if(child.pid()==-sign*413 || child.pid()==-sign*423) {
	  ++nprod;
	  Dstar.push_back(child);
	}
	else if(child.pid()==-sign*15) {
	  ++nprod;
	  tau.push_back(child);
	}
	else if(child.pid()==sign*16) {
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
	  findChildren(child,sign,nprod,Dstar,tau,nu);
	}
      }
    }

    void findTau(const Particle & p, int & sign, unsigned int & nprod,
		 Particles & piP, Particles & pi0, Particles & nu) {
      for(const Particle & child : p.children()) {
	if(child.pid()==111) {
	  ++nprod;
	  pi0.push_back(child);
	}
	else if(child.pid()==sign*211) {
	  ++nprod;
	  piP.push_back(child);
	}
	else if(child.pid()==-sign*16) {
	  ++nprod;
	  nu.push_back(child);
	}
	else if(child.pid()==22)
	  continue;
	else if(child.children().empty() || child.pid()==221 || child.pid()==331) {
	  ++nprod;
	}
	else {
	  findTau(child,sign,nprod,piP,pi0,nu);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over B0 mesons 
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::B0 or
									       Cuts::abspid==PID::BPLUS)) {
	// find the B decay
	int sign = p.pid()/p.abspid();
	unsigned int nprod = 0;
	Particles Dstar,tau,nu;
	findChildren(p,sign,nprod,Dstar,tau,nu);
	if(nprod!=3 || Dstar.size()!=1 || tau.size() !=1 || nu.size()!=1)
	  continue;
	// check decay
	if(p.pid()==PID::B0) {
	  if(Dstar[0].pid()!=-sign*413) vetoEvent;
	}
	else if(p.pid()==PID::BPLUS) {
	  if(Dstar[0].pid()!=-sign*423) vetoEvent;
	}
	// find the tau decay
	nprod=0;
	nu.clear();
	Particles piP,pi0;
	findTau(tau[0],sign,nprod,piP,pi0,nu);
	if(nu.size()!=1) continue;
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	FourMomentum ptau = boost1.transform(tau[0].momentum());
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ptau);
	// pion mode
	if(nprod==2 && piP.size()==1) {
	  FourMomentum pPi = boost2.transform(boost1.transform(piP[0].momentum()));
	  double cTheta = pPi.p3().unit().dot(ptau.p3().unit());
	  _h_pi->fill(cTheta);
	}
	// rho mode
	else if(nprod==3 && piP.size()==1 && pi0.size()==1) {
	  FourMomentum pRho = boost2.transform(boost1.transform(piP[0].momentum()+pi0[0].momentum()));
	  double cTheta = pRho.p3().unit().dot(ptau.p3().unit());
	  _h_rho->fill(cTheta);
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
      normalize(_h_pi );
      normalize(_h_rho);
      // the polarization
      Scatter2DPtr _h_alpha;
      book(_h_alpha,1,1,2);
      pair<double,double> alpha_pi  = calcAlpha(_h_pi );
      pair<double,double> alpha_rho = calcAlpha(_h_rho);
      // 0.45 factor for rho
      alpha_rho.first /=0.46;
      alpha_rho.second/=0.46;
      pair<double,double> alpha;
      alpha.first  = (alpha_pi.first*sqr(alpha_rho.second)+alpha_rho.first*sqr(alpha_pi.second))/(sqr(alpha_pi.second)+sqr(alpha_rho.second));
      alpha.second = alpha_pi.second*alpha_rho.second/sqrt(sqr(alpha_pi.second)+sqr(alpha_rho.second));
      _h_alpha->addPoint(0.5,alpha.first, make_pair(0.5,0.5), make_pair(alpha.second,alpha.second) );

    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_pi,_h_rho;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2018_I1621272);

}
