// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Thrust.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class OPAL_1997_I440103 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_1997_I440103);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(Thrust(FinalState()), "Thrust");
      declare(ChargedFinalState(), "FS");
      declare(InitialQuarks(), "IQF");
      declare(UnstableParticles(), "UFS" );

      // Book histograms
      // B*
      book(_h_B , 8,1,1);
      book(_h_B2, "/TMP/c_theta_B", 20, -1.,1.);
      // phi
      book(_h_phi_ctheta , 5,1,1);
      book(_h_phi_ctheta2, "/TMP/c_theta_phi2", 20, -1.,1.   );
      book(_h_phi_ctheta3, "/TMP/c_theta_phi3", 20, -1.,1.   );
      book(_h_phi_ctheta4, "/TMP/c_theta_phi4", 20, -1.,1.   );
      book(_h_phi_alpha  , 5,1,2);
      book(_h_phi_alpha2 , "/TMP/alpha_phi2", 20, 0.,0.5*M_PI);
      book(_h_phi_alpha3 , "/TMP/alpha_phi3", 20, 0.,0.5*M_PI);
      book(_h_phi_alpha4 , "/TMP/alpha_phi4", 20, 0.,0.5*M_PI);
      book(_h_phi_beta   , 5,1,3);
      book(_h_phi_beta2  , "/TMP/beta_phi2", 20, 0.,0.5*M_PI );
      book(_h_phi_beta3  , "/TMP/beta_phi3", 20, 0.,0.5*M_PI );
      book(_h_phi_beta4  , "/TMP/beta_phi4", 20, 0.,0.5*M_PI );
      book(_c_phi_cos_plus , "/TMP/c_phi_cos_plus1");
      book(_c_phi_cos_neg  , "/TMP/c_phi_cos_neg1" );
      book(_c_phi_sin_plus , "/TMP/c_phi_sin_plus1");
      book(_c_phi_sin_neg  , "/TMP/c_phi_sin_neg1" );
      book(_c_phi_cos_plus2, "/TMP/c_phi_cos_plus2");
      book(_c_phi_cos_neg2 , "/TMP/c_phi_cos_neg2" );
      book(_c_phi_sin_plus2, "/TMP/c_phi_sin_plus2");
      book(_c_phi_sin_neg2 , "/TMP/c_phi_sin_neg2" );
      book(_c_phi_cos_plus3, "/TMP/c_phi_cos_plus3");
      book(_c_phi_cos_neg3 , "/TMP/c_phi_cos_neg3" );
      book(_c_phi_sin_plus3, "/TMP/c_phi_sin_plus3");
      book(_c_phi_sin_neg3 , "/TMP/c_phi_sin_neg3" );
      // D*
      book(_h_DS_ctheta , 6,1,1);
      book(_h_DS_ctheta2, "/TMP/c_theta_DS2", 20, -1.,1.   );
      book(_h_DS_alpha  , 7,1,1);
      book(_h_DS_alpha2 , "/TMP/alpha_DS2", 20, 0.,0.5*M_PI);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // thrust, to define an axis
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      
      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      /// @todo Yuck... does this *really* have to be quark-based?!?
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      } else {
        map<int, double> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }

      // loop over the particles
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==513 or Cuts::abspid==523 or
										Cuts::pid==333    or Cuts::abspid==413)) {
	int sign = p.pid()/p.abspid();
	Particle decay;
	if(p.children().size()!=2) continue;
	// B*
	if(p.abspid()==513 or p.abspid()==523) {
	  int mid = p.abspid()-2;
	  if(p.children()[0].pid()==sign*mid && 
	     p.children()[1].pid()==22) {
	    decay = p.children()[1];
	  }
	  else if(p.children()[1].pid()==sign*mid && 
		  p.children()[0].pid()==22) {
	    decay = p.children()[0];
	  }
	  else
	    continue;
	}
	// phi
	else if(p.pid()==333) {
	  // cut x_E > 0.7
	  double xE = p.momentum().E()/meanBeamMom;
	  if(xE<0.7) continue;
	  if(p.children()[0].pid()== 321 && 
	     p.children()[1].pid()==-321) {
	    decay = p.children()[0];
	  }
	  else if(p.children()[1].pid()== 321 && 
		  p.children()[0].pid()==-321) {
	    decay = p.children()[1];
	  }
	  else
	    continue;
	}
	// D*
	else if(p.abspid()==413) {
	  double xE = p.momentum().E()/meanBeamMom;
	  if(xE<0.5 || flavour!=4) continue;
	  if(p.children()[0].pid()==sign*421 && 
	     p.children()[1].pid()==sign*211) {
	    decay = p.children()[1];
	  }
	  else if(p.children()[1].pid()==sign*421 && 
		  p.children()[0].pid()==sign*211) {
	    decay = p.children()[0];
	  }
	  else
	    continue;
	  
	}
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 e1z = p.p3().unit();	
	FourMomentum pp = boost.transform(decay.momentum());
	Vector3 axis1 = boost.transform(decay.momentum()).p3().unit();
	double ctheta = e1z.dot(axis1);
	if(p.abspid()==513 or p.abspid()==523) {
	  _h_B ->fill(ctheta);
	  _h_B2->fill(ctheta);
	}
	// D*
	else if(p.abspid()==413) {
	  // y and z axis
	  Vector3 e1y = e1z.cross(axis).unit();
	  Vector3 e1x = e1y.cross(e1z).unit();
	  // helicity beam axis, all phis
	  // cos theta_H
	  _h_DS_ctheta ->fill(ctheta);
	  _h_DS_ctheta2->fill(ctheta);
	  // alpha
	  double phi = atan2(e1y.dot(axis1),e1x.dot(axis1));
	  double alpha = abs(abs(phi)-0.5*M_PI);
	  _h_DS_alpha   ->fill(alpha);
	  _h_DS_alpha2  ->fill(alpha);
	}
	else if(p.pid()==333) {
	  // y and z axis
	  Vector3 e1y = e1z.cross(axis).unit();
	  Vector3 e1x = e1y.cross(e1z).unit();
	  // helicity beam axis, all phis
	  // cos theta_H
	  _h_phi_ctheta->fill(abs(ctheta));
	  _h_phi_ctheta2->fill(ctheta);
	  // alpha and beta
	  double phi = atan2(e1y.dot(axis1),e1x.dot(axis1));
	  double alpha = abs(abs(phi)-0.5*M_PI);
	  double beta =  abs(abs(phi+0.25*M_PI)-0.5*M_PI);
	  _h_phi_alpha   ->fill(alpha);
	  _h_phi_alpha2  ->fill(alpha);
	  _h_phi_beta    ->fill( beta);
	  _h_phi_beta2   ->fill( beta);
	  /// counters for asymmetries
	  double sin2H = 2.*ctheta*sqrt(1.-sqr(ctheta));
	  if(sin2H*cos(phi)>0.) 
	    _c_phi_cos_plus->fill();
	  else
	    _c_phi_cos_neg->fill();
	  if(sin2H*sin(phi)>0.) 
	    _c_phi_sin_plus->fill();
	  else
	    _c_phi_sin_neg->fill();
	  // whether or not is a primary hadron
	  Particle parent = p.parents()[0];
	  if(parent.children().size()==1 && parent.abspid()==p.abspid())
	    parent = parent.parents()[0];
	  bool primary = !PID::isHadron(parent.pid());
	  if(primary) {
	    // cos theta_H
	    _h_phi_ctheta3->fill(ctheta);
	    // alpha and beta
	    _h_phi_alpha3  ->fill(alpha);
	    _h_phi_beta3   ->fill( beta);
	    /// counters for asymmetries
	    if(sin2H*cos(phi)>0.) 
	      _c_phi_cos_plus2->fill();
	    else
	      _c_phi_cos_neg2->fill();
	    if(sin2H*sin(phi)>0.) 
	      _c_phi_sin_plus2->fill();
	    else
	      _c_phi_sin_neg2->fill();
	  }
	  // pT w.r.t thrust axis
	  double pT = sqrt(sqr(thrust.thrustMajorAxis().dot(p.momentum().p3()))+
			   sqr(thrust.thrustMinorAxis().dot(p.momentum().p3())));
	  // helicity-quark frame
	  if(pT>1.2) {
	    // cos theta H
	    _h_phi_ctheta4->fill(ctheta);
	    Vector3 axis2;
	    if(p.momentum().p3().dot(thrust.thrustAxis())>=0.) {
	      axis2 = thrust.thrustAxis();
	    }
	    else {
	      axis2 =-thrust.thrustAxis();
	    }
	    Vector3 e2y = e1z.cross(axis2).unit();
	    Vector3 e2x = e2y.cross(e1z).unit();
	    // alpha and beta
	    double phi = atan2(e2y.dot(axis1),e2x.dot(axis1));
	    double alpha = abs(abs(phi)-0.5*M_PI);
	    double beta =  abs(abs(phi+0.25*M_PI)-0.5*M_PI);
	    _h_phi_alpha4  ->fill(alpha);
	    _h_phi_beta4   ->fill( beta);
	    /// counters for asymmetries
	    double sin2H = 2.*ctheta*sqrt(1.-sqr(ctheta));
	    if(sin2H*cos(phi)>0.) 
	      _c_phi_cos_plus3->fill();
	    else
	      _c_phi_cos_neg3->fill();
	    if(sin2H*sin(phi)>0.) 
	      _c_phi_sin_plus3->fill();
	    else
	      _c_phi_sin_neg3->fill();
	    
	  }
	}
      }
    }
    
    pair<double,double> calcRho(Histo1DPtr hist,unsigned int mode) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
      double sum1(0.),sum2(0.);
      for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai(0.),bi(0.);
	if(mode==0) {
	  ai = 0.25*( -bin.xMin()*(3.-sqr(bin.xMin())) + bin.xMax()*(3.-sqr(bin.xMax())));
	  bi =-0.75*( -bin.xMin()*(1.-sqr(bin.xMin())) + bin.xMax()*(1.-sqr(bin.xMax())));
	}
	else if(mode==1) {
	  ai = 0.125*( -bin.xMin()*(3.+sqr(bin.xMin())) + bin.xMax()*(3.+sqr(bin.xMax())));
	  bi = 0.375*( -bin.xMin()*(1.-sqr(bin.xMin())) + bin.xMax()*(1.-sqr(bin.xMax())));
	}
	else if(mode==2) {
	  ai = -2.*(bin.xMin()-bin.xMax())/M_PI;
	  bi = -2.*(sin(2.*bin.xMin())-sin(2.*bin.xMax()))/M_PI;
	}
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // B*
      normalize(_h_B,1.,false);
      normalize(_h_B2);
      pair<double,double> rho = calcRho(_h_B2,1);
      Scatter2DPtr h_rhoB;
      book(h_rhoB,4,1,1);
      h_rhoB->addPoint(91.2, rho.first, make_pair(0.5,0.5),
		       make_pair(rho.second,rho.second) );
      // D*
      normalize(_h_DS_ctheta );
      normalize(_h_DS_ctheta2);
      rho = calcRho(_h_DS_ctheta2,1);
      Scatter2DPtr h_rhoD;
      book(h_rhoD,3,1,1);
      h_rhoD->addPoint(1., rho.first, make_pair(0.5,0.5),
			make_pair(rho.second,rho.second) );
      normalize(_h_DS_alpha );
      normalize(_h_DS_alpha2);
      Scatter2DPtr h_reRho_D;
      book(h_reRho_D,3,1,2);
      rho = calcRho(_h_DS_alpha2,2);
      h_reRho_D->addPoint(1., rho.first, make_pair(0.5,0.5),
			  make_pair(rho.second,rho.second) );
      // phi
      // rho00
      normalize(_h_phi_ctheta );
      normalize(_h_phi_ctheta2);
      normalize(_h_phi_ctheta3);
      normalize(_h_phi_ctheta4);
      Scatter2DPtr hrho_phi;
      book(hrho_phi,1,1,1);
      rho = calcRho(_h_phi_ctheta2,0);
      hrho_phi->addPoint(1., rho.first, make_pair(0.5,0.5),
			 make_pair(rho.second,rho.second) );
      rho = calcRho(_h_phi_ctheta3,0);
      hrho_phi->addPoint(2., rho.first, make_pair(0.5,0.5),
			 make_pair(rho.second,rho.second) );
      rho = calcRho(_h_phi_ctheta4,0);
      hrho_phi->addPoint(3., rho.first, make_pair(0.5,0.5),
			 make_pair(rho.second,rho.second) );
      // Re rho
      normalize(_h_phi_alpha );
      normalize(_h_phi_alpha2);
      normalize(_h_phi_alpha3);
      normalize(_h_phi_alpha4);
      Scatter2DPtr  hreRho_phi;
      book(hreRho_phi,1,1,2);
      rho = calcRho(_h_phi_alpha2,2);
      hreRho_phi->addPoint(1., rho.first, make_pair(0.5,0.5),
			   make_pair(rho.second,rho.second) );
      rho = calcRho(_h_phi_alpha3,2);
      hreRho_phi->addPoint(2., rho.first, make_pair(0.5,0.5),
			   make_pair(rho.second,rho.second) );
      rho = calcRho(_h_phi_alpha4,2);
      hreRho_phi->addPoint(3., rho.first, make_pair(0.5,0.5),
			   make_pair(rho.second,rho.second) );
      // Im rho
      normalize(_h_phi_beta );
      normalize(_h_phi_beta2);
      normalize(_h_phi_beta3);
      normalize(_h_phi_beta4);
      Scatter2DPtr himRho_phi;
      book(himRho_phi,1,1,3);
      rho = calcRho(_h_phi_beta2,2);
      himRho_phi->addPoint(1., rho.first, make_pair(0.5,0.5),
			   make_pair(rho.second,rho.second) );
      rho = calcRho(_h_phi_beta3,2);
      himRho_phi->addPoint(2., rho.first, make_pair(0.5,0.5),
			   make_pair(rho.second,rho.second) );
      rho = calcRho(_h_phi_beta4,2);
      himRho_phi->addPoint(3., rho.first, make_pair(0.5,0.5),
			   make_pair(rho.second,rho.second) );
      // real diff
      Scatter1D temp = (*_c_phi_cos_plus-*_c_phi_cos_neg)/(*_c_phi_cos_plus+*_c_phi_cos_neg);
      Scatter1D temp2 = (*_c_phi_cos_plus2-*_c_phi_cos_neg2)/(*_c_phi_cos_plus2+*_c_phi_cos_neg2);
      Scatter1D temp3 = (*_c_phi_cos_plus3-*_c_phi_cos_neg3)/(*_c_phi_cos_plus3+*_c_phi_cos_neg3);
      Scatter2DPtr hreDiff_phi;
      book(hreDiff_phi,1,1,4);
      hreDiff_phi->addPoint(1., temp.points()[0].x(), make_pair(0.5,0.5),
			    make_pair(temp.points()[0].xErrMinus(),temp.points()[0].xErrPlus()) );
      hreDiff_phi->addPoint(2., temp2.points()[0].x(), make_pair(0.5,0.5),
			    make_pair(temp2.points()[0].xErrMinus(),temp2.points()[0].xErrPlus()) );
      hreDiff_phi->addPoint(3., temp3.points()[0].x(), make_pair(0.5,0.5),
			    make_pair(temp3.points()[0].xErrMinus(),temp3.points()[0].xErrPlus()) );
      // im diff
      temp  = (*_c_phi_sin_plus-*_c_phi_sin_neg)/(*_c_phi_sin_plus+*_c_phi_sin_neg);
      temp2 = (*_c_phi_sin_plus2-*_c_phi_sin_neg2)/(*_c_phi_sin_plus2+*_c_phi_sin_neg2);
      temp3 = (*_c_phi_sin_plus3-*_c_phi_sin_neg3)/(*_c_phi_sin_plus3+*_c_phi_sin_neg3);
      Scatter2DPtr himDiff_phi;
      book(himDiff_phi,1,1,5);
      himDiff_phi->addPoint(1., temp.points()[0].x(), make_pair(0.5,0.5),
			    make_pair(temp.points()[0].xErrMinus(),temp.points()[0].xErrPlus()) );
      himDiff_phi->addPoint(2., temp2.points()[0].x(), make_pair(0.5,0.5),
			    make_pair(temp2.points()[0].xErrMinus(),temp2.points()[0].xErrPlus()) );
      himDiff_phi->addPoint(3., temp3.points()[0].x(), make_pair(0.5,0.5),
			    make_pair(temp3.points()[0].xErrMinus(),temp3.points()[0].xErrPlus()) );
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_B,_h_B2;
    Histo1DPtr _h_phi_ctheta, _h_phi_ctheta2, _h_phi_ctheta3, _h_phi_ctheta4;
    Histo1DPtr _h_phi_alpha , _h_phi_alpha2 , _h_phi_alpha3 , _h_phi_alpha4 ;
    Histo1DPtr _h_phi_beta  , _h_phi_beta2  , _h_phi_beta3  , _h_phi_beta4  ;
    CounterPtr _c_phi_cos_plus, _c_phi_cos_neg, _c_phi_cos_plus2, _c_phi_cos_neg2, _c_phi_cos_plus3, _c_phi_cos_neg3;
    CounterPtr _c_phi_sin_plus, _c_phi_sin_neg, _c_phi_sin_plus2, _c_phi_sin_neg2, _c_phi_sin_plus3, _c_phi_sin_neg3;
    Histo1DPtr _h_DS_ctheta, _h_DS_ctheta2;
    Histo1DPtr _h_DS_alpha , _h_DS_alpha2 ;
    //@}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OPAL_1997_I440103);


}
