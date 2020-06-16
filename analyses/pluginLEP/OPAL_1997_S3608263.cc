// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief OPAL K*0 fragmentation function paper
  /// @author Peter Richardson
  class OPAL_1997_S3608263 : public Analysis {
  public:

    /// Constructor
    OPAL_1997_S3608263()
      : Analysis("OPAL_1997_S3608263")
    {}


    /// @name Analysis methods
    //@{
    
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_histXeK0   , 1, 1, 1);
      {Histo1DPtr temp; _h_ctheta.add(0.   ,0.01 ,book(temp, "ctheta_00",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.01 ,0.03 ,book(temp, "ctheta_01",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.03 ,0.10 ,book(temp, "ctheta_02",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.10 ,0.125,book(temp, "ctheta_03",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.125,0.14 ,book(temp, "ctheta_04",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.14 ,0.16 ,book(temp, "ctheta_05",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.16 ,0.20 ,book(temp, "ctheta_06",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.20 ,0.30 ,book(temp, "ctheta_07",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.30 ,0.40 ,book(temp, "ctheta_08",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.40 ,0.50 ,book(temp, "ctheta_09",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.50 ,0.70 ,book(temp, "ctheta_10",20,-1.,1.));}
      {Histo1DPtr temp; _h_ctheta.add(0.70 ,1.00 ,book(temp, "ctheta_11",20,-1.,1.));}

      book(_h_ctheta_large,"ctheta_large",20,-1.,1.);

      {Histo1DPtr temp; _h_alpha.add(0.   ,0.01 ,book(temp, "alpha_00",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.01 ,0.03 ,book(temp, "alpha_01",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.03 ,0.10 ,book(temp, "alpha_02",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.10 ,0.125,book(temp, "alpha_03",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.125,0.14 ,book(temp, "alpha_04",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.14 ,0.16 ,book(temp, "alpha_05",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.16 ,0.20 ,book(temp, "alpha_06",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.20 ,0.30 ,book(temp, "alpha_07",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.30 ,0.40 ,book(temp, "alpha_08",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.40 ,0.50 ,book(temp, "alpha_09",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.50 ,0.70 ,book(temp, "alpha_10",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha.add(0.70 ,1.00 ,book(temp, "alpha_11",20,0.,0.5*M_PI));}

      {Histo1DPtr temp; _h_alpha_low.add(0.30 ,0.40 ,book(temp, "alpha_low_00",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha_low.add(0.40 ,0.50 ,book(temp, "alpha_low_01",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha_low.add(0.50 ,0.70 ,book(temp, "alpha_low_02",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha_low.add(0.70 ,1.00 ,book(temp, "alpha_low_03",20,0.,0.5*M_PI));}
     
      {Histo1DPtr temp; _h_alpha_high.add(0.30 ,0.40 ,book(temp, "alpha_high_00",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha_high.add(0.40 ,0.50 ,book(temp, "alpha_high_01",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha_high.add(0.50 ,0.70 ,book(temp, "alpha_high_02",20,0.,0.5*M_PI));}
      {Histo1DPtr temp; _h_alpha_high.add(0.70 ,1.00 ,book(temp, "alpha_high_03",20,0.,0.5*M_PI));}
     

      book(_h_alpha_large     ,"_h_alpha_large"     ,20,0.,0.5*M_PI);
      book(_h_alpha_large_low ,"_h_alpha_large_low" ,20,0.,0.5*M_PI);
      book(_h_alpha_large_high,"_h_alpha_large_high",20,0.,0.5*M_PI);
    }

    pair<double,double> calcRho(Histo1DPtr hist,unsigned int imode) {
      if(hist->numEntries()==0.) return make_pair(0.,0.);
     double sum1(0.),sum2(0.);
     for (auto bin : hist->bins() ) {
	double Oi = bin.area();
	if(Oi==0.) continue;
	double ai,bi;
	if(imode==0) {
	  ai = 0.25*(bin.xMax()*(3.-sqr(bin.xMax())) - bin.xMin()*(3.-sqr(bin.xMin())));
	  bi = 0.75*(bin.xMin()*(1.-sqr(bin.xMin())) - bin.xMax()*(1.-sqr(bin.xMax())));
	}
	else {
	  ai = 2.*(bin.xMax()-bin.xMin())/M_PI;
	  bi = 2.*(sin(2.*bin.xMax())-sin(2.*bin.xMin()))/M_PI;
	}
	double Ei = bin.areaErr();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
     }
     return make_pair(sum2/sum1,sqrt(1./sum1));
   }
    
    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableFinalState>(e, "UFS");

      for (const Particle& p : ufs.particles(Cuts::abspid==313)) {
	double xp = p.p3().mod()/meanBeamMom;
	_histXeK0->fill(xp);
	int sign = p.pid()/313;
	if(p.children().size()!=2) continue;
	Particle kaon;
	if(p.children()[0].pid()==sign*321 && p.children()[1].pid()==-sign*211) {
	  kaon = p.children()[0];
	}
	else if(p.children()[1].pid()==sign*321 && p.children()[0].pid()==-sign*211) {
	  kaon = p.children()[1];
	}
	else
	  continue;
	// spin axes
	Vector3 e1z = p.momentum().p3().unit();
	Vector3 e1y = e1z.cross(axis).unit();
	Vector3 e1x = e1y.cross(e1z).unit();
	LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
	Vector3 axis1 = boost.transform(kaon.momentum()).p3().unit();
	double ctheta = e1z.dot(axis1);
	double phi = atan2(e1y.dot(axis1),e1x.dot(axis1));
	double alpha = abs(abs(phi)-0.5*M_PI);
	_h_ctheta.fill(xp,ctheta);
	_h_alpha .fill(xp,alpha );
	double cBeam = axis.dot(e1z);
	if(cBeam<0.5)
	  _h_alpha_low .fill(xp,alpha);
	else
	  _h_alpha_high.fill(xp,alpha);

	if(xp>0.3) {
	  _h_ctheta_large->fill(ctheta);
	  _h_alpha_large->fill(alpha);
	  if(cBeam<0.5)
	    _h_alpha_large_low ->fill(alpha);
	  else
	    _h_alpha_large_high->fill(alpha);
	}
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXeK0, 1./sumOfWeights());
      vector<double> x = {0., 0.01 ,0.03 ,0.10 ,0.125,0.14 ,0.16 ,0.20 ,0.30 ,0.40 ,0.50 ,0.70 ,1.00};
      Scatter2DPtr h_rho00;
      book(h_rho00,2,1,1);
      Scatter2DPtr h_rho_off;
      book(h_rho_off,2,1,2);
      Scatter2DPtr h_ratio;
      book(h_ratio,3,1,1);
      Scatter2DPtr h_off_low;
      book(h_off_low,4,1,1);
      Scatter2DPtr h_off_high;
      book(h_off_high,4,1,2);
      for(unsigned int ix=0;ix<_h_ctheta.histos().size();++ix) {
	// extract the rho00 component
	normalize(_h_ctheta.histos()[ix]);
	pair<double,double> rho00 = calcRho(_h_ctheta.histos()[ix],0);
	h_rho00->addPoint(0.5*(x[ix]+x[ix+1]), rho00.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			  make_pair(rho00.second,rho00.second) );
	// extract the rho 1 -1 component
	normalize(_h_alpha.histos()[ix]);
	pair<double,double> rho_off = calcRho(_h_alpha.histos()[ix],1);
	h_rho_off->addPoint(0.5*(x[ix]+x[ix+1]), rho_off.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			  make_pair(rho_off.second,rho_off.second) );
	if(ix<8) continue;
	// at large xp also the ratio
	double ratio = rho_off.first/(1.-rho00.first);
	double dr    = ((rho_off.second - rho00.first*rho_off.second + rho_off.first*rho00.second))/sqr(1.-rho00.first);
	h_ratio->addPoint(0.5*(x[ix]+x[ix+1]), ratio, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			  make_pair(dr,dr) );
	unsigned int iy=ix-8;
	// rho 1 -1 for cos theta <0.5
	normalize(_h_alpha_low .histos()[iy]);
	rho_off = calcRho(_h_alpha_low.histos()[iy],1);
	h_off_low->addPoint(0.5*(x[ix]+x[ix+1]), rho_off.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			    make_pair(rho_off.second,rho_off.second) );
	// rho 1 -1 for cos theta >0.5
	normalize(_h_alpha_high.histos()[iy]);
	rho_off = calcRho(_h_alpha_high.histos()[iy],1);
	h_off_high->addPoint(0.5*(x[ix]+x[ix+1]), rho_off.first, make_pair(0.5*(x[ix+1]-x[ix]),0.5*(x[ix+1]-x[ix])),
			    make_pair(rho_off.second,rho_off.second) );
	
      }
      // ratio for xp 0.3
      normalize(_h_ctheta_large);
      pair<double,double> rho00 = calcRho(_h_ctheta_large,0);
      normalize(_h_alpha_large);
      pair<double,double> rho_off = calcRho(_h_alpha_large,1);
      double ratio = rho_off.first/(1.-rho00.first);
      double dr    = ((rho_off.second - rho00.first*rho_off.second + rho_off.first*rho00.second))/sqr(1.-rho00.first);
      Scatter2DPtr h_ratio_large;
      book(h_ratio_large,3,2,1);
      h_ratio_large->addPoint(0.65, ratio, make_pair(0.35,0.35),make_pair(dr,dr) );
      // rho 1 -1 for xp >0.3 and cos theta < 0.5
      normalize(_h_alpha_large_low );
      rho_off = calcRho(_h_alpha_large_low,1);
      book(h_off_low,4,2,1);
      h_off_low->addPoint(0.65, rho_off.first, make_pair(0.35,0.35),make_pair(rho_off.second,rho_off.second) );
      // rho 1 -1 for xp >0.3 and cos theta > 0.5
      normalize(_h_alpha_large_high);
      rho_off = calcRho(_h_alpha_large_high,1);
      book(h_off_high,4,2,2);
      h_off_high->addPoint(0.65, rho_off.first, make_pair(0.35,0.35),make_pair(rho_off.second,rho_off.second) );
    }

    //@}


  private:

    Histo1DPtr _histXeK0;
    BinnedHistogram _h_ctheta,_h_alpha,_h_alpha_low,_h_alpha_high;
    Histo1DPtr  _h_ctheta_large,_h_alpha_large, _h_alpha_large_low,_h_alpha_large_high;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1997_S3608263);

}
