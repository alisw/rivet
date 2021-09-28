// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BESIII_2017_I1506414 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2017_I1506414);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      declare(FinalState(), "FS");
	      
      // Book histograms
      if(fuzzyEquals(sqrtS(),3.1,1e-1)) {
	book(_h_xi , 1, 1, 2);
	book(_h_sig, 1, 1, 1);
      }
      else if (fuzzyEquals(sqrtS(), 3.686, 1E-1)) {
	book(_h_xi , 1, 1, 4);
	book(_h_sig, 1, 1, 3);
      }
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


      const UnstableParticles & ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p :  ufs.particles(Cuts::abspid==3322 or Cuts::abspid==3214)) {
       	if(p.children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(p,nRes,ncount);
	bool matched=false;
	// check for antiparticle
	for (const Particle& p2 :  ufs.particles(Cuts::pid==-p.pid())) {
	  if(p2.children().empty()) continue;
	  map<long,int> nRes2=nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2==0) {
	    matched = true;
	    for(auto const & val : nRes2) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    // fond baryon and antibaryon
	    if(matched) {
	      // calc cosine
	      double ctheta;
	      if(p.pid()>0)
		ctheta = p .momentum().p3().unit().dot(axis);
	      else
		ctheta = p2.momentum().p3().unit().dot(axis);
	      if(abs(p.pid())==3322)
		_h_xi ->fill(ctheta);
	      else if(abs(p.pid())==3214)
		_h_sig->fill(ctheta);
	      break;
	    }
	  }
	}
	if(matched) break;
      }
    }
    
    pair<double,pair<double,double> > calcAlpha(Histo1DPtr hist) {
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
      // find energy
      int ioff=-1;
      if(fuzzyEquals(sqrtS(),3.1,1e-1)) {
	ioff=0;
      }
      else if (fuzzyEquals(sqrtS(), 3.686, 1E-1)) {
	ioff=1;
      }
      normalize(_h_xi,1.,false);
      Scatter2DPtr _h_alpha_xi;
      book(_h_alpha_xi,2,2*ioff+2,1);
      pair<double,pair<double,double> > alpha = calcAlpha(_h_xi);
      _h_alpha_xi->addPoint(0.5, alpha.first, make_pair(0.5,0.5),
			    make_pair(alpha.second.first,alpha.second.second) );
      normalize(_h_sig,1.,false);
      Scatter2DPtr _h_alpha_sig;
      book(_h_alpha_sig,2,2*ioff+1,1);
      alpha = calcAlpha(_h_sig);
      _h_alpha_sig->addPoint(0.5, alpha.first, make_pair(0.5,0.5),
			     make_pair(alpha.second.first,alpha.second.second) );
      
    }
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_xi,_h_sig;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BESIII_2017_I1506414);


}
