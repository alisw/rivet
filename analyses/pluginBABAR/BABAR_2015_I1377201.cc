// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief azimuthal asymmetries in pipi Kpi and KK
  class BABAR_2015_I1377201 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2015_I1377201);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      const FinalState fs;
      declare(fs,"FS");
      declare(Thrust(fs),"Thrust");
      declare(Beam(), "Beams");
      // declare the histos for the distributions
      string type  [3] = {"KK","Kpi","pipi"};
      string charge[3] = {"Like","Opposite","All"};
      unsigned int nbin=20;
      for(unsigned int itype=0;itype<3;++itype) {
	for(unsigned int icharge=0;icharge<3;++icharge) {
	  for(unsigned int ibin=0;ibin<16;++ibin) {
	    std::ostringstream title1;
	    title1 << "/TMP/h_thrust" << type[itype] << "_" << charge[icharge] << "_" << ibin+1;
	    book(_h_thrust[itype][icharge][ibin],title1.str(),nbin,0.,M_PI);
	    std::ostringstream title2;
	    title2 << "/TMP/h_hadron" << type[itype] << "_" << charge[icharge] << "_" << ibin+1;
	    book(_h_hadron[itype][icharge][ibin],title2.str(),nbin,0.,M_PI);
	  }
	}
      }
    }

    unsigned int iBin(double z) {
      if     (z<.2) return 0;
      else if(z<.3) return 1;
      else if(z<.5) return 2;
      else          return 3;
    }
    
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis1;
      if(beams.first.pid()>0)
	axis1 = beams.first .momentum().p3().unit();
      else
	axis1 = beams.second.momentum().p3().unit();
      // apply thrust cuts  T > 0.8  and | cos θ th | < 0.6
      Thrust thrust = apply<Thrust>(event,"Thrust");
      if(thrust.thrust()<=0.8) vetoEvent;
      if(cos(thrust.thrustAxis().polarAngle())>=0.6) vetoEvent;
      // construct x,y,z axes for thrust defn
      ThreeVector t_z = thrust.thrustAxis();
      ThreeVector t_x = (axis1-t_z.dot(axis1)*t_z).unit();
      ThreeVector t_y = t_z.cross(t_x);
      // loop over the particles
      Particles charged = apply<FinalState>(event,"FS").particles(Cuts::abspid==PID::PIPLUS || Cuts::abspid==PID::KPLUS);
      for(unsigned int ix=0;ix<charged.size();++ix) {
	// z and angle cut
	const double x1=2.*charged[ix].momentum().t()/sqrtS();
	if(x1<0.16||x1>.9) continue;
	if(abs(t_z.angle(charged[ix].momentum().p3()))>0.25*M_PI) continue;
	double dot1 = t_z.dot(charged[ix].p3());
	for(unsigned int iy=ix+1;iy<charged.size();++iy) {
	  const double x2=2.*charged[iy].momentum().t()/sqrtS();
	  // z and angle cut
	  if(x2<0.16||x2>.9) continue;
	  if(abs(t_z.angle(charged[ix].momentum().p3()))>0.25*M_PI) continue;
	  // different hemi
	  double dot2 = t_z.dot(charged[iy].p3());
	  if(dot1*dot2>0.) continue;
	  Particle p1=charged[ix], p2=charged[iy];
	  double z1(x1),z2(x2);
	  // randomly order the particles
	  if(rand()/static_cast<double>(RAND_MAX) < 0.5 ) {
	    swap(p1,p2);
	    swap(z1,z2);
	  }
	  // thrust def
	  double phi12 = atan2(p1.p3().dot(t_y),p1.p3().dot(t_x))+atan2(p2.p3().dot(t_y),p2.p3().dot(t_x));
	  if(phi12>M_PI)  phi12 -= 2*M_PI;
	  if(phi12<-M_PI) phi12 += 2*M_PI;
	  if(phi12<0.) phi12 = -phi12;
	  // hadron defn
	  ThreeVector h_z = p2.p3().unit();
	  ThreeVector h_x = (axis1-h_z.dot(axis1)*h_z).unit();
	  ThreeVector pt1 = p1.p3()-h_z.dot(p1.p3())*h_z;
	  double phi0 = pt1.angle(h_x);
	  if(phi0>M_PI)  phi0 -= 2*M_PI;
	  if(phi0<-M_PI) phi0 += 2*M_PI;
	  int ibin = 4*iBin(z1)+iBin(z2);
	  // pi pi
	  if(p1.abspid()==PID::PIPLUS && p2.abspid()==PID::PIPLUS) {
	    if(p1.pid()==p2.pid()) {
	      _h_thrust[2][0][ibin]->fill(phi12);
	      _h_hadron[2][0][ibin]->fill(phi0);
	    }
	    else {
	      _h_thrust[2][1][ibin]->fill(phi12);
	      _h_hadron[2][1][ibin]->fill(phi0);
	    }
	    _h_thrust[2][2][ibin]->fill(phi12);
	    _h_hadron[2][2][ibin]->fill(phi0);
	  }
	  // K K
	  else if(p1.abspid()==PID::KPLUS && p2.abspid()==PID::KPLUS) {
	    if(p1.pid()==p2.pid()) {
	      _h_thrust[0][0][ibin]->fill(phi12);
	      _h_hadron[0][0][ibin]->fill(phi0);
	    }
	    else {
	      _h_thrust[0][1][ibin]->fill(phi12);
	      _h_hadron[0][1][ibin]->fill(phi0);
	    }
	    _h_thrust[0][2][ibin]->fill(phi12);
	    _h_hadron[0][2][ibin]->fill(phi0);
	  }
	  // K pi
	  else {
	    if(p1.pid()*p2.pid()>0) {
	      _h_thrust[1][0][ibin]->fill(phi12);
	      _h_hadron[1][0][ibin]->fill(phi0);
	    }
	    else {
	      _h_thrust[1][1][ibin]->fill(phi12);
	      _h_hadron[1][1][ibin]->fill(phi0);
	    }
	    _h_thrust[1][2][ibin]->fill(phi12);
	    _h_hadron[1][2][ibin]->fill(phi0);
	  }

	  
	}
      }
    }
    
    pair<double,double> calcAsymmetry(Scatter2DPtr hist,double fact=1.) {
      double sum1(0.),sum2(0.);
      for (auto point : hist->points() ) {
	double Oi = point.y();
	if(Oi==0. || std::isnan(Oi) ) continue;
	double ai = 1.;
	double bi = (sin(fact*point.xMax())-sin(fact*point.xMin()))/(point.xMax()-point.xMin())/fact;
	double Ei = point.yErrAvg();
	sum1 += sqr(bi/Ei);
	sum2 += bi/sqr(Ei)*(Oi-ai);
      }
      if(sum1==0.) return make_pair(0.,0.);
      return make_pair(sum2/sum1*1e4,sqrt(1./sum1)*1e4);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int itype=0;itype<3;++itype) {
	for(unsigned int icharge=0;icharge<3;++icharge) {
	  for(unsigned int ibin=0;ibin<16;++ibin) {
	    normalize(_h_thrust[itype][icharge][ibin]);
	    normalize(_h_hadron[itype][icharge][ibin]);
	  }
	}
      }
      // construct ther ratios
      // declare the histos for the distributions
      string type  [3] = {"pipi","Kpi","KK"};
      string charge[3] = {"Like","Opposite","All"};
      for(unsigned int itype=0;itype<3;++itype) {
	Scatter3DPtr h3_thrust_UL;
	book(h3_thrust_UL,2*itype+1,1,2,0);
	Scatter3DPtr h3_thrust_UC;
	book(h3_thrust_UC,2*itype+1,1,3,0);
	Scatter3DPtr h3_hadron_UL;
	book(h3_hadron_UL,2*itype+2,1,2,0);
	Scatter3DPtr h3_hadron_UC;
	book(h3_hadron_UC,2*itype+2,1,3,0);
	
	unsigned int ihist=1;
	Scatter2DPtr h2_thrust_UL;
	book(h2_thrust_UL,7+2*itype,ihist,2);
	Scatter2DPtr h2_thrust_UC;
	book(h2_thrust_UC,7+2*itype,ihist,3);
	Scatter2DPtr h2_hadron_UL;
	book(h2_hadron_UL,8+2*itype,ihist,2);
	Scatter2DPtr h2_hadron_UC;
	book(h2_hadron_UC,8+2*itype,ihist,3);
	
	Scatter3D temphisto1(refData<Scatter3D>(2*itype+1, 1, 2));
	Scatter3D temphisto2(refData<Scatter3D>(2*itype+2, 1, 2));
	for(unsigned int ibin=0;ibin<16;++ibin) {
	  const Point3D & p1 = temphisto1.points()[ibin];
	  const Point3D & p2 = temphisto2.points()[ibin];
	  if(ibin>0 && ibin%4==0) {
	    ++ihist;
	    book(h2_thrust_UL,7+2*itype,ihist,2);
	    book(h2_thrust_UC,7+2*itype,ihist,3);
	    book(h2_hadron_UL,8+2*itype,ihist,2);
	    book(h2_hadron_UC,8+2*itype,ihist,3);
	  }
	  // thrust direction
	  // opposite/like sign
	  std::ostringstream title1;
	  title1 << "/TMP/R_thrust_" << type[itype] << "_UL_" << ibin+1;
	  Scatter2DPtr htemp;
	  book(htemp,title1.str());
	  divide(_h_thrust[itype][1][ibin],
		 _h_thrust[itype][0][ibin],htemp);
	  pair<double,double> asym = calcAsymmetry(htemp);
	  h3_thrust_UL->addPoint(p1.x()    ,p1.y()    ,asym.first,
				 p1.xErrs(),p1.yErrs(),make_pair(asym.second,asym.second) );
	  h2_thrust_UL->addPoint(p1.y()    ,asym.first,p1.yErrs(),make_pair(asym.second,asym.second) );
	  // opposite/all sign
	  std::ostringstream title2;
	  title2 << "/TMP/R_thrust_" << type[itype] << "_UC_" << ibin+1;
	  book(htemp,title2.str());
	  divide(_h_thrust[itype][1][ibin],
		 _h_thrust[itype][2][ibin],htemp);
	  asym = calcAsymmetry(htemp);
	  h3_thrust_UC->addPoint(p1.x()    ,p1.y()    ,asym.first,
				 p1.xErrs(),p1.yErrs(),make_pair(asym.second,asym.second) );
	  h2_thrust_UC->addPoint(p1.y()    ,asym.first,p1.yErrs(),make_pair(asym.second,asym.second) );
	  // hadron dirn
	  // opposite/like sign
	  std::ostringstream title3;
	  title3 << "/TMP/R_hadron_" << type[itype] << "_UL_" << ibin+1;
	  book(htemp,title3.str());
	  divide(_h_hadron[itype][1][ibin],
		 _h_hadron[itype][0][ibin],htemp);
	  asym = calcAsymmetry(htemp,2.);
	  h3_hadron_UL->addPoint(p2.x()    ,p2.y()    ,asym.first,
				 p2.xErrs(),p2.yErrs(),make_pair(asym.second,asym.second) );
	  h2_hadron_UL->addPoint(p2.y()    ,asym.first,p2.yErrs(),make_pair(asym.second,asym.second) );
	  // opposite/all sign
	  std::ostringstream title4;
	  title4 << "/TMP/R_hadron_" << type[itype] << "_UC_" << ibin+1;
	  book(htemp,title4.str());
	  divide(_h_hadron[itype][1][ibin],
		 _h_hadron[itype][2][ibin],htemp);
	  asym = calcAsymmetry(htemp,2.);
	  h3_hadron_UC->addPoint(p2.x()    ,p2.y()    ,asym.first,
				 p2.xErrs(),p2.yErrs(),make_pair(asym.second,asym.second) );
	  h2_hadron_UC->addPoint(p2.y()    ,asym.first,p2.yErrs(),make_pair(asym.second,asym.second) );
	}
      }
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_thrust[3][3][16],_h_hadron[3][3][16];
    //@}


  };


  DECLARE_RIVET_PLUGIN(BABAR_2015_I1377201);

}
