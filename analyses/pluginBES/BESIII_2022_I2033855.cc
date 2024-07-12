// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> gamma chi_c0,2 -> Xi Xibar
  class BESIII_2022_I2033855 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2033855);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(Cuts::pid==10441 || Cuts::pid==445), "UFS");
      declare(FinalState(), "FS");
      // book hists
      for(unsigned int ix=0;ix<3;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],4+ix,1,1+iy);
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.children().empty()) {
	  nRes[child.pid()]-=1;
	  --ncount;
	}
	else
	  findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // get the axis, direction of incoming electron
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      Vector3 axis;
      if(beams.first.pid()>0)
	axis = beams.first .momentum().p3().unit();
      else
	axis = beams.second.momentum().p3().unit();
      // types of final state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // loop over chi_c states
      Particle chi;
      bool matched = false;
      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p :  ufs.particles()) {
       	if(p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
	if(ncount==1) {
	  matched = true;
	  for(auto const & val : nRes) {
	    if(val.first==PID::PHOTON) {
	      if(val.second!=1) {
	      matched = false;
	      break;
	      }
	    }
	    else if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    chi=p;
	    break;
	  }
	}
      }
      if(!matched) vetoEvent;
      // have chi_c find psi2S 
      if(chi.parents().empty() || chi.children().size()!=2 ||
	 chi.children()[0].pid() != -chi.children()[1].pid()) vetoEvent;
      Particle psi2S = chi.parents()[0];
      if(psi2S.pid()!=100443 || psi2S.children().size()!=2) vetoEvent;
      // then the first photon
      Particle gamma1;
      if(psi2S.children()[0].pid()==PID::PHOTON)
	gamma1 = psi2S.children()[0];
      else if(psi2S.children()[1].pid()==PID::PHOTON)
	gamma1 = psi2S.children()[1];
      else
	vetoEvent;
      // now the decay products of the chi_c
      Particle bPlus,bMinus;
      bool foundBaryon=false;
      for(unsigned int ix=0;ix<2;++ix) {
	if(chi.children()[ix].pid()==PID::XIMINUS ||
	   chi.children()[ix].pid()==PID::XI0 ) {
	  foundBaryon=true;
	  bPlus=chi.children()[ix];
	}
	else if(chi.children()[ix].pid()==-PID::XIMINUS ||
		chi.children()[ix].pid()==-PID::XI0 ) {
	  bMinus=chi.children()[ix];
	}
      }
      if(!foundBaryon) vetoEvent;
      // type chi state
      unsigned int ichi = 0;
      if(chi.pid()==20443) ichi = 1;
      else if(chi.pid()==445) ichi = 2;
      LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(chi.momentum().betaVec());
      Vector3 e1z = gamma1.momentum().p3().unit();
      FourMomentum pBaryon = boost1.transform(bPlus.momentum());
      Vector3 axis1 = pBaryon.p3().unit();
      double cBaryon = e1z.dot(axis1);
      if(bPlus.pid()==PID::XIMINUS)
	_h[ichi][0]->fill(cBaryon);
      else
	_h[ichi][1]->fill(cBaryon);
    }
    
    pair<double,pair<double,double> > calcAlpha0(Histo1DPtr hist) {
      if(hist->numEntries()==0.) return make_pair(0.,make_pair(0.,0.));
      double d = 3./(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double c = 3.*(hist->xMax()-hist->xMin())/(pow(hist->xMax(),3)-pow(hist->xMin(),3));
      double sum1(0.),sum2(0.),sum3(0.),sum4(0.),sum5(0.);
      for (auto bin : hist->bins() ) {
       	double Oi = bin.area();
	if(Oi==0.) continue;
	double a =  d*(bin.xMax() - bin.xMin());
	double b = d/3.*(pow(bin.xMax(),3) - pow(bin.xMin(),3));
       	double Ei = bin.areaErr();
	sum1 +=   a*Oi/sqr(Ei);
	sum2 +=   b*Oi/sqr(Ei);
	sum3 += sqr(a)/sqr(Ei);
	sum4 += sqr(b)/sqr(Ei);
	sum5 +=    a*b/sqr(Ei);
      }
      // calculate alpha
      double alpha = (-c*sum1 + sqr(c)*sum2 + sum3 - c*sum5)/(sum1 - c*sum2 + c*sum4 - sum5);
      // and error
      double cc = -pow((sum3 + sqr(c)*sum4 - 2*c*sum5),3);
      double bb = -2*sqr(sum3 + sqr(c)*sum4 - 2*c*sum5)*(sum1 - c*sum2 + c*sum4 - sum5);
      double aa =  sqr(sum1 - c*sum2 + c*sum4 - sum5)*(-sum3 - sqr(c)*sum4 + sqr(sum1 - c*sum2 + c*sum4 - sum5) + 2*c*sum5);      
      double dis = sqr(bb)-4.*aa*cc;
      if(dis>0.) {
	dis = sqrt(dis);
	return make_pair(alpha,make_pair(0.5*(-bb+dis)/aa,-0.5*(-bb-dis)/aa));
      }
      else {
	return make_pair(alpha,make_pair(0.,0.));
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<3;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  normalize(_h[ix][iy],1.,false);
	  pair<double,pair<double,double> > alpha0 = calcAlpha0(_h[ix][iy]);
	  Scatter2DPtr _h_alpha0;
	  book(_h_alpha0,1+ix,1,1+iy);
	  _h_alpha0->addPoint(0.5, alpha0.first, make_pair(0.5,0.5),
			      make_pair(alpha0.second.first,alpha0.second.second) );
	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2033855);

}
