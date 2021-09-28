// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include <sstream>
namespace Rivet {


  /// @brief Lambda polarization at LEP1
  class OPAL_1997_I447188 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_1997_I447188);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      const ChargedFinalState cfs;
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(UnstableParticles(), "UFS");

      // Book the histograms
      {Histo1DPtr temp; _h_ctheta.add(0.027, 0.05 ,book(temp, "/TMP/ctheta_0",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.05 , 0.08 ,book(temp, 4,1,1));}
      {Histo1DPtr temp; _h_ctheta.add(0.08 , 0.09 ,book(temp, "/TMP/ctheta_2",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.09 , 0.1  ,book(temp, "/TMP/ctheta_3",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.1  , 0.15 ,book(temp, 4,1,2));}
      {Histo1DPtr temp; _h_ctheta.add(0.15 , 0.2  ,book(temp, "/TMP/ctheta_5",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.2  , 0.3  ,book(temp, 4,1,3));}
      {Histo1DPtr temp; _h_ctheta.add(0.3  , 0.4  ,book(temp, "/TMP/ctheta_7",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.4  , 1.0  ,book(temp, "/TMP/ctheta_8",20,-1.,1.));}
      book(_h_ctheta_large, 4,1,4);
      
      {Histo1DPtr temp; _h_cphi .add(0.3,0.6,book(temp, "/TMP/cphiP_0",10,0.,1.));}
      {Histo1DPtr temp; _h_cphi .add(0.3,0.6,book(temp, 5,1,1));}
      {Histo1DPtr temp; _h_cphi .add(0.6,0.9,book(temp, 5,1,2));}
      {Histo1DPtr temp; _h_cphi .add(0.9,1.2,book(temp, "/TMP/cphiP_3",10,0.,1.));}
      {Histo1DPtr temp; _h_cphi .add(1.2,1.5,book(temp, "/TMP/cphiP_4",10,0.,1.));}
      
      book(_h_cphi_low , 5,1,4);
      book(_h_cphi_mid , "/TMP/cphiP_mid" ,10,0.,1.);
      book(_h_cphi_high, 5,1,3);

      {Histo1DPtr temp; _h_plus_lam.add(0.027,0.05,book(temp, "/TMP/lamP_0",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.05 ,0.08,book(temp, "/TMP/lamP_1",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.08 ,0.09,book(temp, "/TMP/lamP_2",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.09 ,0.1 ,book(temp, "/TMP/lamP_3",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.1  ,0.15,book(temp, "/TMP/lamP_4",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.15 ,0.2 ,book(temp, "/TMP/lamP_5",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.2  ,0.3 ,book(temp, "/TMP/lamP_6",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.3  ,0.4 ,book(temp, "/TMP/lamP_7",20,-1.,1.));}
      {Histo1DPtr temp; _h_plus_lam.add(0.4  ,1.0 ,book(temp, "/TMP/lamP_8",20,-1.,1.));}
      
      {Histo1DPtr temp; _h_minus_lam.add(0.027,0.05,book(temp, "/TMP/lamM_0",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.05 ,0.08,book(temp, "/TMP/lamM_1",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.08 ,0.09,book(temp, "/TMP/lamM_2",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.09 ,0.1 ,book(temp, "/TMP/lamM_3",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.1  ,0.15,book(temp, "/TMP/lamM_4",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.15 ,0.2 ,book(temp, "/TMP/lamM_5",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.2  ,0.3 ,book(temp, "/TMP/lamM_6",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.3  ,0.4 ,book(temp, "/TMP/lamM_7",20,-1.,1.));}
      {Histo1DPtr temp; _h_minus_lam.add(0.4  ,1.0 ,book(temp, "/TMP/lamM_8",20,-1.,1.));}

      {Histo1DPtr temp; _h_plus_lam_large1  = book(temp, "/TMP/lamP_large_1",20,-1.,1.);}
      {Histo1DPtr temp; _h_plus_lam_large2  = book(temp, "/TMP/lamP_large_2",20,-1.,1.);}
      {Histo1DPtr temp; _h_minus_lam_large1 = book(temp, "/TMP/lamM_large_1",20,-1.,1.);}
      {Histo1DPtr temp; _h_minus_lam_large2 = book(temp, "/TMP/lamM_large_2",20,-1.,1.);}
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      Vector3 beamAxis;
      if(beams.first.pid()==-11) {
	beamAxis = beams.first .momentum().p3().unit();
      }
      else {
	beamAxis = beams.second.momentum().p3().unit();
      }
	
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // thrust, to define an axis
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const UnstableParticles& ufs = apply<UnstableFinalState>(event, "UFS");

      for(const Particle & lambda : ufs.particles(Cuts::abspid==3122)) {
	double xE = lambda.momentum().t()/meanBeamMom;
	int sign = lambda.pid()/3122;
	Vector3 axis1 = lambda.momentum().p3().unit();
	// assymetry
	double cLam = axis1.dot(beamAxis);
	if(sign>0)
	  _h_plus_lam .fill(xE,cLam);
	else
	  _h_minus_lam.fill(xE,cLam);
	if(xE>0.15) {
	  if(sign>0)
	    _h_plus_lam_large1 ->fill(cLam);
	  else
	    _h_minus_lam_large1->fill(cLam);
	}
	if(xE>0.3) {
	  if(sign>0)
	    _h_plus_lam_large2 ->fill(cLam);
	  else
	    _h_minus_lam_large2->fill(cLam);
	}
	if(lambda.children().size()!=2) continue;
	// look at the decay products
	Particle proton,pion;
	if(lambda.children()[0].pid()==sign*2212 && 
	   lambda.children()[1].pid()==-sign*211) {
	  proton = lambda.children()[0];
	  pion   = lambda.children()[1];
	}
	else if(lambda.children()[1].pid()==sign*2212 && 
		lambda.children()[0].pid()==-sign*211) {
	  proton = lambda.children()[1];
	  pion   = lambda.children()[0];
	}
	else
	  continue;
	// boost to the lambda rest frame
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(lambda.momentum().betaVec());
	FourMomentum pproton = boost.transform(proton.momentum());
	// longitudinal polarization
	double ctheta = axis1.dot(pproton.p3().unit());
	_h_ctheta.fill(xE,ctheta);
	if(xE>=0.3)  _h_ctheta_large->fill(ctheta);
	// transverse polarization
	Vector3 axis2;
	if(lambda.momentum().p3().dot(thrust.thrustAxis())>=0.) {
	  axis2 = thrust.thrustAxis();
	}
	else {
	  axis2 =-thrust.thrustAxis();
	}
	Vector3 axis3 = axis2.cross(axis1).unit();	
	double pT = sqrt(sqr(thrust.thrustMajorAxis().dot(lambda.momentum().p3()))+
			 sqr(thrust.thrustMinorAxis().dot(lambda.momentum().p3())));
	double cPhi = axis3.dot(pproton.p3().unit());
	_h_cphi.fill(pT,cPhi);
	if(pT>0.3) _h_cphi_low ->fill(cPhi);
	if(pT>0.6) _h_cphi_mid ->fill(cPhi);
	if(pT>1.5) _h_cphi_high->fill(cPhi);
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
    
    pair<double,double> calcAsymmetry(Scatter2DPtr hist,unsigned int mode) {
      double sum1(0.),sum2(0.);
      for (auto bin : hist->points() ) {
	double Oi = bin.y();
	if(Oi==0.) continue;
	double bi;
	if(mode==0)
	  bi = 0.25*(bin.xMax()-bin.xMin())*(bin.xMax()+bin.xMin());
	else
	  bi = 4.*(bin.xMax()+bin.xMin())/(3.+sqr(bin.xMax())+bin.xMax()*bin.xMin()+sqr(bin.xMin()));
	double Ei = bin.yErrAvg();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*Oi;
      }
      return make_pair(sum2/sum1,sqrt(1./sum1));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // longitudinal polarization
      vector<double> x_val = {3.850000e-02,6.500000e-02,8.500000e-02,9.500000e-02,1.250000e-01,
			      1.750000e-01,2.500000e-01,3.500000e-01,7.000000e-01};
      vector<double> x_err = {1.150000e-02,1.500000e-02,5.000000e-03,5.000000e-03,2.500000e-02,
			      2.500000e-02,5.000000e-02,5.000000e-02,3.000000e-01};
      unsigned int ipoint=0;
      double aLam = 0.642;
      Scatter2DPtr h_long;
      book(h_long, 1,1,1);
      for(Histo1DPtr & hist : _h_ctheta.histos()) {
	normalize(hist);
	pair<double,double> alpha = calcAlpha(hist);
	alpha.first  /=aLam;
	alpha.second /=aLam;
	h_long->addPoint(x_val[ipoint], 100.*alpha.first, make_pair(x_err[ipoint],x_err[ipoint]),
			   make_pair(100.*alpha.second,100.*alpha.second) );	
	++ipoint;
      }
      normalize(_h_ctheta_large);
      pair<double,double> alpha = calcAlpha(_h_ctheta_large);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Scatter2DPtr h_long_l;
      book(h_long_l,1,1,2);
      h_long_l->addPoint(0.65, 100.*alpha.first, make_pair(0.35,0.35), make_pair(100.*alpha.second,100.*alpha.second) );
      // transverse polarization
      vector<double> pT_val = {0.15,0.45,0.75,1.05,1.35};
      vector<double> pT_err = {0.15,0.15,0.15,0.15,0.15};
      Scatter2DPtr h_trans;
      book(h_trans,2,1,1);
      ipoint=0;
      for(Histo1DPtr & hist : _h_cphi.histos()) {
	normalize(hist);
	pair<double,double> alpha = calcAlpha(hist);
	alpha.first  /=aLam;
	alpha.second /=aLam;
	h_trans->addPoint(pT_val[ipoint], 100.*alpha.first, make_pair(pT_err[ipoint],pT_err[ipoint]),
			   make_pair(100.*alpha.second,100.*alpha.second) );	
	++ipoint;
      }
      normalize(_h_cphi_low);
      alpha = calcAlpha(_h_cphi_low);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Scatter2DPtr h_trans_low;
      book(h_trans_low,2,1,2);
      h_trans_low->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(100.*alpha.second,100.*alpha.second) );
      normalize(_h_cphi_mid);
      alpha = calcAlpha(_h_cphi_mid);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Scatter2DPtr h_trans_mid;
      book(h_trans_mid,2,1,3);
      h_trans_mid->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(100.*alpha.second,100.*alpha.second) );
      normalize(_h_cphi_high);
      alpha = calcAlpha(_h_cphi_high);
      alpha.first  /=aLam;
      alpha.second /=aLam;
      Scatter2DPtr h_trans_high;
      book(h_trans_high,2,1,4);
      h_trans_high->addPoint(0.5, alpha.first, make_pair(0.5,0.5), make_pair(100.*alpha.second,100.*alpha.second) );
      // asyymetry
      Scatter2DPtr h_asym;
      book(h_asym,3,1,1);
      for(unsigned int ix=0;ix<9;++ix) {
       	normalize(_h_plus_lam.histos()[ix] );
       	normalize(_h_minus_lam.histos()[ix]);
	std::ostringstream title;
       	title << "/TMP/a_lam_" << ix;
       	Scatter2DPtr sTemp;
	book(sTemp,title.str());
       	asymm(_h_plus_lam.histos()[ix],_h_minus_lam.histos()[ix],sTemp);
       	pair<double,double> alpha = calcAsymmetry(sTemp,1);
       	h_asym->addPoint(x_val[ix],-alpha.first, make_pair(x_err[ix],x_err[ix]),
			 make_pair(alpha.second,alpha.second) );
      }
      normalize(_h_plus_lam_large1 );
      normalize(_h_minus_lam_large1);
      Scatter2DPtr sTemp;
      book(sTemp,"/TMP/a_lam_large1");
      asymm(_h_plus_lam_large1,_h_minus_lam_large1,sTemp);
      alpha = calcAsymmetry(sTemp,1);
      book(h_asym,3,1,2);
      h_asym->addPoint(0.575,-alpha.first, make_pair(0.425,0.425),
		       make_pair(alpha.second,alpha.second) );
      normalize(_h_plus_lam_large2 );
      normalize(_h_minus_lam_large2);
      book(sTemp,"/TMP/a_lam_large2");
      asymm(_h_plus_lam_large2,_h_minus_lam_large2,sTemp);
      alpha = calcAsymmetry(sTemp,1);
      book(h_asym,3,1,3);
      h_asym->addPoint(0.65,-alpha.first, make_pair(0.35,0.35),
		       make_pair(alpha.second,alpha.second) );
    }

    //@}


    /// @name Histograms
    //@{
    BinnedHistogram _h_ctheta,_h_cphi,_h_plus_lam,_h_minus_lam;
    Histo1DPtr _h_ctheta_large;
    Histo1DPtr _h_cphi_low, _h_cphi_mid, _h_cphi_high;
    Histo1DPtr _h_plus_lam_large1,_h_plus_lam_large2;
    Histo1DPtr _h_minus_lam_large1,_h_minus_lam_large2;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1997_I447188);


}
