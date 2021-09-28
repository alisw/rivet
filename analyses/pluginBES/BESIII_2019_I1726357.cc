// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda Lambdabar cross section
  class BESIII_2019_I1726357 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BESIII_2019_I1726357);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_sigma ,1,1,1);
      book(_h_cTheta,2,1,1);
      double xlow=-1., step=0.2;
      for(unsigned int ix=0;ix<10;++ix) {
	Histo1DPtr temp;
	std::ostringstream title;
	title << "/TMP/h_pol_" << ix;
	book(temp,title.str(),20,-1.,1.);
	_h_pol.add(xlow,xlow+step,temp);
	xlow+=step;
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // find the Lambdas
      bool matched = false;
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      Particle Lambda;
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
	const Particle& p1 = ufs.particles()[ix];
	if(abs(p1.pid())!=3122) continue;
	// check fs
	bool fs = true;
	for (const Particle & child : p1.children()) {
	  if(child.pid()==p1.pid()) {
	    fs = false;
	    break;
	  }
	}
	if(!fs) continue;
	// find the children
	map<long,int> nRes = nCount;
	int ncount = ntotal;
	findChildren(p1,nRes,ncount);
	for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
	  matched=false;
	  const Particle& p2 = ufs.particles()[iy];
	  if(abs(p2.pid())!=3122) continue;
	  // check fs
	  bool fs = true;
	  for (const Particle & child : p2.children()) {
	    if(child.pid()==p2.pid()) {
	      fs = false;
	      break;
	    }
	  }
	  if(!fs) continue;
	  map<long,int> nRes2 = nRes;
	  int ncount2 = ncount;
	  findChildren(p2,nRes2,ncount2);
	  if(ncount2!=0) continue;
	  matched=true;
	  for(auto const & val : nRes2) {
	    if(val.second!=0) {
	      matched = false;
	      break;
	    }
	  }
	  if(matched) {
	    _h_sigma->fill(2.396);
	    if(p1.pid()==PID::LAMBDA) {
	      Lambda=p1;
	    }
	    else {
	      Lambda=p2;
	    }
	    break;
	  }
	}
	if(matched) break;
      }
      // now for the polarization measurements
      if(matched) {
	double cTheta = Lambda.momentum().p3().unit().dot(axis);
	_h_cTheta->fill(cTheta);
	Particle proton;
	if(Lambda.children().size()!=2) return;
	if(Lambda.children()[0].pid()==PID::PROTON &&
	   Lambda.children()[1].pid()==PID::PIMINUS)
	  proton = Lambda.children()[0];
	else if(Lambda.children()[1].pid()==PID::PROTON &&
		Lambda.children()[0].pid()==PID::PIMINUS)
	  proton = Lambda.children()[1];
	else return;
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(Lambda.momentum().betaVec());
	Vector3 axis1 = boost1.transform(proton.momentum()).p3().unit();
	double cPhi = axis1.dot(Lambda.momentum().p3().unit());
	_h_pol.fill(cTheta,cPhi);
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
      scale(_h_sigma,crossSection()/sumOfWeights()/picobarn);
      normalize(_h_cTheta);
      Scatter2DPtr _h_alpha;
      book(_h_alpha,3,1,1);
      double step=0.2, x=-0.9;
      for(unsigned int ix=0;ix<10;++ix) {
	normalize(_h_pol.histos()[ix]);
	pair<double,double> alpha = calcAlpha(_h_pol.histos()[ix]);
	_h_alpha->addPoint(x, alpha.first, make_pair(0.5*step,0.5*step), make_pair(alpha.second,alpha.second) );
	x+=step;
      }
    }
    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_sigma,_h_cTheta;
    BinnedHistogram _h_pol;
    ///@}


  };


  DECLARE_RIVET_PLUGIN(BESIII_2019_I1726357);

}
