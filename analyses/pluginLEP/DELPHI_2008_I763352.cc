// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief tau polarization at 197 GeV
  class DELPHI_2008_I763352 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(DELPHI_2008_I763352);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      book(_h_e  ,"_t_e " , 20,-1,1);
      book(_h_mu ,"_t_mu" , 20,-1,1);
      book(_h_pi ,"_t_pi" , 20,-1,1);
      book(_h_rho,"_t_rho", 20,-1,1);
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
      // loop over tau leptons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==15)) {
	unsigned int nprod(0);
	Particles piP, pi0, ell, nu_ell, nu_tau;
	findTau(p,nprod,piP, pi0, ell, nu_ell, nu_tau);
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	if(nprod==2 && nu_tau.size()==1 && piP.size()==1) {
	  FourMomentum pPi = boost1.transform(piP[0].momentum());
	  double cTheta = pPi.p3().unit().dot(p.momentum().p3().unit());
	  _h_pi->fill(cTheta);
	}
	else if(nprod==3 && nu_tau.size()==1 && ell.size()==1 && nu_ell.size()==1) {
	  if(ell[0].pid()==PID::ELECTRON) {
	    _h_e ->fill(2.*ell[0].momentum().t()/sqrtS());
	  }
	  else {
	    _h_mu->fill(2.*ell[0].momentum().t()/sqrtS());
	  }
	}
	else if(nprod==3 && nu_tau.size()==1 && piP.size()==1&& pi0.size()==1) {
	  FourMomentum pRho = boost1.transform(piP[0].momentum()+pi0[0].momentum());
	  double cTheta = pRho.p3().unit().dot(p.momentum().p3().unit());
	  _h_rho->fill(cTheta);
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
      Scatter2DPtr h_P;
      book(h_P,1,1,1);
      normalize(_h_e);
      pair<double,double> P_e  = calcP(_h_e,1);
      double s1 = P_e.first/sqr(P_e.second);
      double s2 = 1./sqr(P_e.second);
      normalize(_h_mu);
      pair<double,double> P_mu = calcP(_h_mu,1);
      s1 += P_mu.first/sqr(P_mu.second);
      s2 += 1./sqr(P_mu.second);
      normalize(_h_pi);
      pair<double,double> P_pi = calcP(_h_pi,0);
      s1 += P_pi.first/sqr(P_pi.second);
      s2 += 1./sqr(P_pi.second);
      normalize(_h_rho);
      pair<double,double> P_rho = calcP(_h_rho,0);
      s1 += P_rho.first/sqr(P_rho.second);
      s2 += 1./sqr(P_rho.second);
      P_rho.first  /=0.46;
      P_rho.second /=0.46;
      // average
      pair<double,double> P_aver = make_pair(s1/s2,sqrt(1./s2));
      h_P->addPoint(0.5,P_aver.first, make_pair(0.5,0.5),
		    make_pair(P_aver.second,P_aver.second));
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_e,_h_mu,_h_pi,_h_rho;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(DELPHI_2008_I763352);

}
