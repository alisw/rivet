// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+e- > tau+ tau-
  class OPAL_2001_I554583 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2001_I554583);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      double xmin=-0.9;
      double step=0.18;
      for(unsigned int ix=0;ix<10;++ix) {
	Histo1DPtr temp;
	std::ostringstream title1;
	title1 << "_h_e_" << ix;
	book(temp,title1.str(), 20,-1,1);
	_h_e .add(xmin,xmin+step,temp);
	std::ostringstream title2;
	title2 << "_h_mu_" << ix;
	book(temp,title2.str(), 20,-1,1);
	_h_mu.add(xmin,xmin+step,temp);
	std::ostringstream title3;
	title3 << "_h_pi_" << ix;
	book(temp,title3.str(), 20,-1,1);
	_h_pi.add(xmin,xmin+step,temp);
	std::ostringstream title4;
	title4 << "_h_rho_" << ix;
	book(temp,title4.str(), 20,-1,1);
	_h_rho.add(xmin,xmin+step,temp);
	xmin+=step;
      }
      book(_t_e  ,"_t_e " , 20,-1,1);
      book(_t_mu ,"_t_mu" , 20,-1,1);
      book(_t_pi ,"_t_pi" , 20,-1,1);
      book(_t_rho,"_t_rho", 20,-1,1);
    }

    void findTau(const Particle & p, unsigned int & nprod,
     		 Particles & piP,Particles & pi0, Particles & ell, Particles & nu_ell,
		 Particles & nu_tau) {
      for(const Particle & child : p.children()) {
	if(child.pid()==PID::ELECTRON || child.pid()==PID::MUON) {
	  ++nprod;
	  ell.push_back(child);
	}
	else if(child.pid()==PID::NU_EBAR || child.pid()==PID::NU_MUBAR) {
	  ++nprod;
	  nu_ell.push_back(child);
	}
	else if(child.pid()==PID::PIMINUS) {
	  ++nprod;
	  piP.push_back(child);
	}
	else if(child.pid()==PID::PI0) {
	  ++nprod;
	  pi0.push_back(child);
	}
	else if(child.pid()==PID::NU_TAU) {
	  ++nprod;
	  nu_tau.push_back(child);
	}
	else if(child.pid()==PID::GAMMA)
	  continue;
	else if(child.children().empty() || child.pid()==221 || child.pid()==331) {
	  ++nprod;
	}
	else {
	  findTau(child,nprod,piP,pi0,ell,nu_ell,nu_tau);
	}
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // require 2 chanrged particles to veto hadronic events
      if(apply<ChargedFinalState>(event, "FS").particles().size()!=2) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // loop over tau leptons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==15)) {
	unsigned int nprod(0);
	Particles piP, pi0, ell, nu_ell, nu_tau;
	findTau(p,nprod,piP, pi0, ell, nu_ell, nu_tau);
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	double cBeam = axis.dot(p.momentum().p3().unit());
	if(nprod==2 && nu_tau.size()==1 && piP.size()==1) {
	  FourMomentum pPi = boost1.transform(piP[0].momentum());
	  double cTheta = pPi.p3().unit().dot(p.momentum().p3().unit());
	  _h_pi. fill(cBeam,cTheta);
	  _t_pi->fill(cTheta);
	}
	else if(nprod==3 && nu_tau.size()==1 && ell.size()==1 && nu_ell.size()==1) {
	  if(ell[0].pid()==PID::ELECTRON) {
	    _h_e . fill(cBeam,2.*ell[0].momentum().t()/sqrtS());
	    _t_e ->fill(2.*ell[0].momentum().t()/sqrtS());
	  }
	  else {
	    _h_mu. fill(cBeam,2.*ell[0].momentum().t()/sqrtS());
	    _t_mu->fill(2.*ell[0].momentum().t()/sqrtS());
	  }
	}
	else if(nprod==3 && nu_tau.size()==1 && piP.size()==1&& pi0.size()==1) {
	  FourMomentum pRho = boost1.transform(piP[0].momentum()+pi0[0].momentum());
	  double cTheta = pRho.p3().unit().dot(p.momentum().p3().unit());
	  _h_rho. fill(cBeam,cTheta);
	  _t_rho->fill(cTheta);
	}
      }
    }

    pair<double,double> calcP(Histo1DPtr hist,unsigned int imode) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai(0.),bi(0.);
	// tau -> pi/rho nu
	if(imode==0) {
	  ai = 0.5*(bin.xMax()-bin.xMin());
	  bi = 0.5*ai*(bin.xMax()+bin.xMin());
	}
	// lepton mode
	else {
	  ai = (-5*bin.xMin() + 3*pow(bin.xMin(),3) -   pow(bin.xMin(),4) + 5*bin.xMax() - 3*pow(bin.xMax(),3) +   pow(bin.xMax(),4))/3.;
	  bi = (  -bin.xMin() + 3*pow(bin.xMin(),3) - 2*pow(bin.xMin(),4) +   bin.xMax() - 3*pow(bin.xMax(),3) + 2*pow(bin.xMax(),4))/3.;
	}
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      Scatter2DPtr _h_P;
      book(_h_P,2,1,1);
      Scatter2DPtr _t_P;
      book(_t_P,1,1,5);
      double x    =-0.81;
      double step = 0.18;
      for(unsigned int ix=0;ix<11;++ix) {
	Histo1DPtr he = ix<10 ? _h_e  .histos()[ix] : _t_e;
      	normalize(he);
       	pair<double,double> P_e  = calcP(he,1);
       	double s1 = P_e.first/sqr(P_e.second);
       	double s2 = 1./sqr(P_e.second);
	Histo1DPtr hmu = ix<10 ? _h_mu  .histos()[ix] : _t_mu;
      	normalize(hmu);
      	pair<double,double> P_mu = calcP(hmu,1);
      	s1 += P_mu.first/sqr(P_mu.second);
      	s2 += 1./sqr(P_mu.second);	
	Histo1DPtr hpi = ix<10 ? _h_pi  .histos()[ix] : _t_pi;
      	normalize(hpi);
  	pair<double,double> P_pi = calcP(hpi,0);
  	s1 += P_pi.first/sqr(P_pi.second);
  	s2 += 1./sqr(P_pi.second);
	Histo1DPtr hrho = ix<10 ? _h_rho  .histos()[ix] : _t_rho;
      	normalize(hrho);
      	pair<double,double> P_rho = calcP(hrho,0);
      	s1 += P_rho.first/sqr(P_rho.second);
      	s2 += 1./sqr(P_rho.second);
      	P_rho.first  /=0.46;
      	P_rho.second /=0.46;
       	// average
      	pair<double,double> P_aver = make_pair(s1/s2,sqrt(1./s2));
	if(ix<10)
	  _h_P->addPoint(x,P_aver.first, make_pair(0.5*step,0.5*step),
			 make_pair(P_aver.second,P_aver.second));
	else
	  _t_P->addPoint(91.2,P_aver.first, make_pair(0.5,0.5),
			 make_pair(P_aver.second,P_aver.second));
	x+=step;
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistogram _h_e,_h_mu,_h_pi,_h_rho;
    Histo1DPtr      _t_e,_t_mu,_t_pi,_t_rho;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(OPAL_2001_I554583);

}
